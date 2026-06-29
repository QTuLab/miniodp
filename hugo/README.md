# Hugo Static Site

`hugo/` builds the public landing page, help page, and species information
pages. It does not run as a long-lived service; build the static files and let
Nginx or another web server serve the generated `public/` directory.

## Main Files

- `config/hugo_default.toml`: default site parameters and route prefixes.
- `data/species_display.toml`: species cards, default links, and public
  statistics.
- `content/species/*.md`: species information pages.
- `static/images/`: species images and shared static assets.

## Build

Use the included default config and set the deployment URL with `--baseURL`.

```bash
cd hugo
hugo --minify --config config/hugo_default.toml --baseURL https://example.org/miniodp/
```

For local preview:

```bash
cd hugo
hugo server --config config/hugo_default.toml --baseURL http://localhost:1313/miniodp/ --port 1313
```

The repository does not include private environment override files. If your
deployment needs different hostnames, analytics, or branding, keep those
override files outside the public repository or pass values through Hugo flags.

## Adding a Species

1. Add the species image under `static/images/`.
2. Add one table to `data/species_display.toml`.
3. Add a matching page under `content/species/<species_key>.md`.
4. Build locally and check the species card, the species page, and all external
   links.

The `species_key` used in `species_display.toml` must match the Dash route and
the species content filename.

## Statistics

The home page displays four public numbers:

- `bulk_datasets`: bulk sequencing runs.
- `bulk_bases`: sequenced bases from bulk assays.
- `single_cell_datasets`: single-cell sample groups.
- `cells`: cells in published single-cell datasets.

Track counts and detailed run-level audit tables are useful for internal
maintenance, but they are not shown on the home page by default.

