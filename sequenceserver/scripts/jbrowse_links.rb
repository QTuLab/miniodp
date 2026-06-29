require 'csv'
require 'json'

module SequenceServer
  module Links
    MINIODP_BASE_PATH = ENV.fetch('MINIODP_BASE_PATH', '/miniodp').sub(%r{/\z}, '').freeze
    MINIODP_JBROWSE_PATH = ('/' + ENV.fetch('MINIODP_JBROWSE_PATH', '/jbrowse2').sub(%r{\A/+}, '').sub(%r{/+\z}, '')).freeze
    MINIODP_BLAST_JBROWSE_PADDING = ENV.fetch('MINIODP_BLAST_JBROWSE_PADDING', '100').to_i
    MINIODP_BLAST_JBROWSE_LINKS = ENV.fetch(
      'MINIODP_BLAST_JBROWSE_LINKS',
      File.expand_path('../config/jbrowse_genome_links.csv', __dir__)
    ).freeze

    JBROWSE_LINKS_BY_FILE = {}
    JBROWSE_LINKS_BY_TITLE = {}

    if File.file?(MINIODP_BLAST_JBROWSE_LINKS)
      CSV.foreach(MINIODP_BLAST_JBROWSE_LINKS, headers: true) do |row|
        enabled = row['enabled'].to_s.strip.downcase
        next unless enabled.empty? || %w[1 true yes y].include?(enabled)

        record = {
          species_key: row['species_key'].to_s.strip,
          assembly: row['assembly'].to_s.strip,
          annotation_tracks: row['annotation_tracks'].to_s.split(',').map(&:strip).reject(&:empty?),
          refname_strip_prefixes: row['refname_strip_prefixes'].to_s.split(',').map(&:strip).reject(&:empty?)
        }
        next if record[:assembly].empty?

        file_name = row['file_name'].to_s.strip
        database_title = row['database_title'].to_s.strip
        JBROWSE_LINKS_BY_FILE[file_name] = record unless file_name.empty?
        JBROWSE_LINKS_BY_TITLE[database_title] = record unless database_title.empty?
      end
    end

    def jbrowse
      record = miniodp_jbrowse_link_record
      return nil unless record

      refname = miniodp_jbrowse_refname(record)
      coords = miniodp_jbrowse_hit_coords
      return nil if refname.empty? || coords.nil?

      start_pos, end_pos = miniodp_jbrowse_view_range(coords)
      locus = "#{refname}:#{start_pos}..#{end_pos}"
      params = [
        ['assembly', record[:assembly]],
        ['loc', locus]
      ]
      tracks = (miniodp_jbrowse_annotation_tracks(record) + [miniodp_jbrowse_blast_track_id]).uniq
      params << ['tracks', tracks.join(',')] unless tracks.empty?
      params << ['sessionTracks', miniodp_jbrowse_session_tracks(refname)]
      url = "#{MINIODP_BASE_PATH}#{MINIODP_JBROWSE_PATH}/?#{miniodp_jbrowse_query(params)}"

      {
        order: 1,
        title: 'JBrowse',
        url: url,
        icon: 'fa-external-link'
      }
    end

    private

    def miniodp_jbrowse_link_record
      records = miniodp_jbrowse_query_database_records
      return records.first if records.length == 1

      records = miniodp_jbrowse_hit_database_records
      return nil unless records.length == 1

      records.first
    rescue StandardError
      nil
    end

    def miniodp_jbrowse_record_for_database(database)
      file_name = File.basename(database.name.to_s)
      title = database.title.to_s
      JBROWSE_LINKS_BY_FILE[file_name] || JBROWSE_LINKS_BY_TITLE[title]
    end

    def miniodp_jbrowse_query_database_records
      records = report.querydb.map do |database|
        miniodp_jbrowse_record_for_database(database)
      end.compact
      records.uniq { |record| [record[:species_key], record[:assembly]] }
    rescue StandardError
      []
    end

    def miniodp_jbrowse_hit_database_records
      records = []
      miniodp_jbrowse_hit_databases.each do |database|
        record = miniodp_jbrowse_record_for_database(database)
        records << record if record
      end
      records.uniq { |record| [record[:species_key], record[:assembly]] }
    end

    def miniodp_jbrowse_hit_databases
      databases = []
      miniodp_jbrowse_database_id_candidates.each do |candidate|
        SequenceServer::Database.each do |database|
          begin
            databases << database if database.include?(candidate)
          rescue StandardError
            next
          end
        end
      end
      databases.uniq { |database| database.name }
    end

    def miniodp_jbrowse_database_id_candidates
      [id, accession, title.to_s.split.first].compact.flat_map do |value|
        miniodp_jbrowse_refname_candidates(value, nil)
      end.uniq
    end

    def miniodp_jbrowse_hit_id_candidates(record)
      values = [title, id, accession, title.to_s.split.first].compact
      values.flat_map do |value|
        miniodp_jbrowse_refname_candidates(value, record)
      end.uniq
    end

    def miniodp_jbrowse_refname_candidates(value, record)
      raw = value.to_s.strip
      return [] if raw.empty?

      candidates = []

      raw.scan(/(?:chromosome|scaffold|contig|primary_assembly):[^:\s]+:([^:\s]+):\d+:\d+:[+-]?\d+/i) do |match|
        candidates << match[0]
      end
      raw.scan(/\b(?:refseq_seqid|seqid|chromosome|scaffold|contig)=([^\s;,]+)/i) do |match|
        candidates << match[0]
      end

      first_token = raw.split.first.to_s
      if first_token.match?(/\A(?:ref|gb|emb|dbj|lcl)\|[^|]+\|?\z/i)
        candidates << first_token.split('|')[1]
        candidates << first_token
        return candidates.compact.uniq
      end

      miniodp_jbrowse_stripped_refnames(first_token, record).each do |candidate|
        candidates << candidate
      end
      candidates << first_token unless first_token.empty?
      if first_token.include?('|')
        first_token.split('|').each do |part|
          part = part.strip
          next if part.empty?
          next if part.match?(/\A(?:ref|gb|emb|dbj|lcl|gnl|BL_ORD_ID)\z/i)

          miniodp_jbrowse_stripped_refnames(part, record).each do |candidate|
            candidates << candidate
          end
          candidates << part
        end
      end
      candidates.uniq
    end

    def miniodp_jbrowse_hit_coords
      points = hsps.flat_map { |hsp| [hsp.sstart, hsp.send] }.compact
      return nil if points.empty?

      [points.min, points.max]
    end

    def miniodp_jbrowse_session_tracks(refname)
      [
        {
          type: 'FeatureTrack',
          trackId: miniodp_jbrowse_blast_track_id,
          name: 'BLAST HSPs',
          assemblyNames: [miniodp_jbrowse_link_record[:assembly]],
          adapter: {
            type: 'FromConfigAdapter',
            features: miniodp_jbrowse_hsp_features(refname)
          },
          displays: [
            {
              type: 'LinearBasicDisplay',
              displayId: "#{miniodp_jbrowse_blast_track_id}-LinearBasicDisplay"
            }
          ]
        }
      ].to_json
    end

    def miniodp_jbrowse_annotation_tracks(record)
      record[:annotation_tracks].to_a
    end

    def miniodp_jbrowse_hsp_features(refname)
      hsps.map.with_index(1) do |hsp, index|
        hsp_start = [hsp.sstart, hsp.send].min
        hsp_end = [hsp.sstart, hsp.send].max
        {
          uniqueId: "blast_hsp_#{index}",
          refName: refname,
          start: hsp_start - 1,
          end: hsp_end,
          name: "BLAST HSP #{index}",
          type: 'match',
          strand: hsp.sstart <= hsp.send ? 1 : -1
        }
      end
    end

    def miniodp_jbrowse_blast_track_id
      'blast_hsps'
    end

    def miniodp_jbrowse_view_range(coords)
      hit_start, hit_end = coords
      hit_span = [hit_end - hit_start + 1, 1].max
      flank = [(hit_span * 0.5).ceil, MINIODP_BLAST_JBROWSE_PADDING].max
      [[hit_start - flank, 1].max, hit_end + flank]
    end

    def miniodp_jbrowse_refname(record)
      miniodp_jbrowse_hit_id_candidates(record).each do |value|
        miniodp_jbrowse_refname_candidates(value, record).each do |candidate|
          next if candidate.empty?
          next if candidate.match?(/\Agnl\|BL_ORD_ID\|\d+\z/)
          next if miniodp_jbrowse_technical_header_token?(candidate)

          return candidate
        end
      end
      ''
    end

    def miniodp_jbrowse_stripped_refnames(value, record)
      return [] unless record

      prefixes = record[:refname_strip_prefixes].to_a.sort_by { |prefix| -prefix.length }
      prefixes.each_with_object([]) do |prefix, stripped|
        stripped << value.delete_prefix(prefix) if value.start_with?(prefix)
      end
    end

    def miniodp_jbrowse_technical_header_token?(value)
      value.match?(/\A(?:dna|dna_rm|dna_sm):(?:chromosome|scaffold|contig|primary_assembly)\z/i)
    end

    def miniodp_jbrowse_query(params)
      params.map do |key, value|
        "#{encode(key)}=#{encode(value)}"
      end.join('&')
    end
  end
end
