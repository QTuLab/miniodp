#!/usr/bin/env ruby

require 'erb'
require_relative 'jbrowse_links'

module SequenceServer
  class Database
    @databases = []

    class << self
      attr_accessor :databases

      def each(&block)
        databases.each(&block)
      end
    end
  end
end

DatabaseStub = Struct.new(:name, :title, :ids) do
  def include?(sequence_id)
    ids.include?(sequence_id)
  end
end

HspStub = Struct.new(:sstart, :send)

HitStub = Struct.new(:id, :accession, :title, :hsps) do
  include SequenceServer::Links

  def encode(value)
    ERB::Util.url_encode(value)
  end
end

def assert(condition, message)
  raise "FAIL: #{message}" unless condition
end

def assert_genome_link(species_label:, database_file:, database_title:, hit_id:, hit_title:, hsp_start:, hsp_end:, expected_assembly:, expected_refname:, expected_track:)
  SequenceServer::Database.databases = [
    DatabaseStub.new(
      "/db/#{database_file}",
      database_title,
      [hit_id]
    )
  ]
  link = HitStub.new(hit_id, hit_id, hit_title, [HspStub.new(hsp_start, hsp_end)]).jbrowse
  assert(link, "#{species_label} JBrowse link")
  assert(link[:url].include?("assembly=#{expected_assembly}"), "#{species_label} assembly")
  assert(link[:url].include?("loc=#{expected_refname}%3A"), "#{species_label} refname")
  assert(!link[:url].include?('loc=dna_rm%3A'), "#{species_label} skip dna_rm token")
  assert(!link[:url].include?('loc=dna%3A'), "#{species_label} skip dna token")
  assert(link[:url].include?(expected_track), "#{species_label} annotation track")
  assert(link[:url].include?('sessionTracks='), "#{species_label} session track")
  link
end

SequenceServer::Database.databases = [
  DatabaseStub.new(
    '/db/Danio_rerio.GRCz11.dna_rm.primary_assembly.fa',
    'Zebrafish Genome, GRCz11',
    ['14']
  )
]
link = HitStub.new('lcl|14', '14', '14 zebrafish chromosome', [HspStub.new(200, 100)]).jbrowse
assert(link[:title] == 'JBrowse', 'JBrowse title')
assert(link[:url].start_with?('/miniodp/jbrowse2/?'), 'root-relative URL')
assert(link[:url].include?('assembly=Danio_rerio'), 'assembly mapping')
assert(link[:url].include?('loc=14%3A1..300'), 'reverse coordinates and 200 percent view range')
assert(link[:url].include?('tracks=Dr_Ensembl100.gff3%2Cblast_hsps'), 'annotation and BLAST HSP tracks')
assert(link[:url].include?('sessionTracks='), 'session track parameter')
assert(
  HitStub.new('lcl|14', '14', '14', [HspStub.new(1, 10)]).send(
    :miniodp_jbrowse_annotation_tracks,
    { assembly: 'Danio_rerio', annotation_tracks: [] }
  ).empty?,
  'annotation tracks come only from CSV config'
)

SequenceServer::Database.databases = [
  DatabaseStub.new(
    '/db/Danio_rerio.GRCz11.dna_rm.primary_assembly.fa',
    'Zebrafish Genome, GRCz11',
    ['Danio_rerio-23', 'dna_rm:primary_assembly']
  )
]
link = HitStub.new(
  'Danio_rerio-23',
  'Danio_rerio-23',
  'Danio_rerio-23 dna_rm:chromosome chromosome:GRCz11:23:1:46223584:1 REF',
  [HspStub.new(7_380_845, 7_381_204)]
).jbrowse
assert(link[:url].include?('assembly=Danio_rerio'), 'zebrafish title assembly mapping')
assert(link[:url].include?('loc=23%3A7380665..7381384'), 'zebrafish title chromosome refname and 200 percent view range')
assert(link[:url].include?('%22refName%22%3A%2223%22'), 'session feature refname')

link = HitStub.new(
  'Danio_rerio-23',
  'Danio_rerio-23',
  'Danio_rerio-23',
  [HspStub.new(7_380_607, 7_381_204)]
).jbrowse
assert(link[:url].include?('loc=23%3A7380308..7381503'), 'zebrafish compact chromosome refname')

link = HitStub.new(
  'dna_rm:primary_assembly',
  'dna_rm:primary_assembly',
  'dna_rm:primary_assembly primary_assembly:GRCz11:23:1:46223584:1 REF',
  [HspStub.new(7_380_607, 7_381_204)]
).jbrowse
assert(link[:url].include?('loc=23%3A7380308..7381503'), 'primary_assembly header refname')
assert(!link[:url].include?('loc=dna_rm%3Aprimary_assembly'), 'skip technical header token')

assert_genome_link(
  species_label: 'cattle',
  database_file: 'Bos_taurus.ARS-UCD1.3.dna_rm.toplevel.fa',
  database_title: 'Cattle Genome, ARS-UCD1.3',
  hit_id: 'Bos_taurus-1',
  hit_title: 'Bos_taurus-1 dna_rm:chromosome chromosome:ARS-UCD1.3:1:1:158534110:1 REF',
  hsp_start: 22_155_470,
  hsp_end: 22_172_856,
  expected_assembly: 'Bos_taurus',
  expected_refname: '1',
  expected_track: 'Bt_Ensembl.gff3'
)

assert_genome_link(
  species_label: 'hydra',
  database_file: 'Hydra_vulgaris_105_v3.dna_rm.toplevel.fa',
  database_title: 'Hydra Genome, v3',
  hit_id: 'Hydra_vulgaris_105_v3-scaffold_1',
  hit_title: 'Hydra_vulgaris_105_v3-scaffold_1 dna_rm:scaffold scaffold:Hydra_105_v3:scaffold_1:1:1000000:1 REF',
  hsp_start: 1_000,
  hsp_end: 1_500,
  expected_assembly: 'Hydra_vulgaris',
  expected_refname: 'scaffold_1',
  expected_track: 'Hv_Ensembl.gff3'
)

assert_genome_link(
  species_label: 'lancelet',
  database_file: 'Branchiostoma_lanceolatum.BraLan2.dna_rm.toplevel.fa',
  database_title: 'Lancelet Genome, BraLan2',
  hit_id: 'Branchiostoma_lanceolatum-BraLan2_1',
  hit_title: 'Branchiostoma_lanceolatum-BraLan2_1 dna_rm:scaffold scaffold:BraLan2:BraLan2_1:1:1000000:1 REF',
  hsp_start: 2_000,
  hsp_end: 2_600,
  expected_assembly: 'Branchiostoma_lanceolatum',
  expected_refname: 'BraLan2_1',
  expected_track: 'Bl_Ensembl.gff3'
)

assert_genome_link(
  species_label: 'medaka',
  database_file: 'medaka_ens94plus.fa',
  database_title: 'Medaka Genome, Ensembl 94',
  hit_id: 'dna_rm:primary_assembly',
  hit_title: 'dna_rm:primary_assembly primary_assembly:ASM223467v1:1:1:37713152:1 REF',
  hsp_start: 22_155_470,
  hsp_end: 22_172_856,
  expected_assembly: 'Oryzias_latipes',
  expected_refname: '1',
  expected_track: 'Ol_Ensembl94.gff3'
)

assert_genome_link(
  species_label: 'mexican tetra',
  database_file: 'Astyanax_mexicanus-2.0.dna_rm.toplevel.fa',
  database_title: 'Mexican tetra Genome, 2.0',
  hit_id: 'Astyanax_mexicanus-2.0-1',
  hit_title: 'Astyanax_mexicanus-2.0-1 dna_rm:chromosome chromosome:Astyanax_mexicanus-2.0:1:1:1000000:1 REF',
  hsp_start: 3_000,
  hsp_end: 3_800,
  expected_assembly: 'Astyanax_mexicanus',
  expected_refname: '1',
  expected_track: 'Am_Ensembl.gff3'
)

assert_genome_link(
  species_label: 'killifish',
  database_file: 'Nothobranchius_furzeri.NfurGRZ-RIMD1.genome.fa',
  database_title: 'Turquoise killifish Genome, NfurGRZ-RIMD1',
  hit_id: 'NC_091749.1',
  hit_title: 'NC_091749.1 Nothobranchius furzeri strain GRZ-AD chromosome 9',
  hsp_start: 26_585_592,
  hsp_end: 26_598_806,
  expected_assembly: 'Nothobranchius_furzeri',
  expected_refname: 'NC_091749.1',
  expected_track: 'Nf_RefSeq.gff3'
)

assert_genome_link(
  species_label: 'killifish RefSeq attribute',
  database_file: 'Nothobranchius_furzeri.NfurGRZ-RIMD1.genome.fa',
  database_title: 'Turquoise killifish Genome, NfurGRZ-RIMD1',
  hit_id: 'refseq_seqid=NC_091755.1',
  hit_title: 'refseq_seqid=NC_091755.1 Nothobranchius furzeri strain GRZ-AD chromosome 15',
  hsp_start: 44_003_471,
  hsp_end: 44_024_589,
  expected_assembly: 'Nothobranchius_furzeri',
  expected_refname: 'NC_091755.1',
  expected_track: 'Nf_RefSeq.gff3'
)

SequenceServer::Database.databases = [
  DatabaseStub.new('/db/Danio_rerio.GRCz11.cdna.all.fa', 'Zebrafish cDNA, Ensembl 115 (GRCz11)', ['ENSDART000001'])
]
assert(HitStub.new('ENSDART000001', 'ENSDART000001', 'ENSDART000001', [HspStub.new(1, 50)]).jbrowse.nil?, 'no cDNA link')

SequenceServer::Database.databases = [
  DatabaseStub.new(
    '/db/Danio_rerio.GRCz11.dna_rm.primary_assembly.fa',
    'Zebrafish Genome, GRCz11',
    ['14']
  ),
  DatabaseStub.new(
    '/db/Nothobranchius_furzeri.NfurGRZ-RIMD1.genome.fa',
    'Turquoise killifish Genome, NfurGRZ-RIMD1',
    ['14']
  )
]
assert(HitStub.new('14', '14', '14 duplicated id', [HspStub.new(10, 20)]).jbrowse.nil?, 'ambiguous genome link')

puts 'JBrowse link tests passed.'
