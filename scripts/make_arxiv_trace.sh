#!/usr/bin/env bash
# Build the arXiv upload of the paper on trace characterisations, man/trace_envelopes/arxiv_trace_inequalities_contrast.zip (and arxiv_trace_inequalities_contrast.tar.gz):
# trace_inequalities_contrast.tex, the generated numbers and the figures (from scripts/trace_characterizations.py and scripts/bounds_search.py in
# man/article, figures/gardner_contrasts.pdf from scripts/trace_envelopes_figures.py, figures/mc_*.pdf from
# scripts/trace_montecarlo.py).
# The archive is compiled twice in a clean directory first; the script fails if pdflatex reports an error or an
# undefined reference. It also writes the compiled man/trace_envelopes/trace_inequalities_contrast.pdf.
set -euo pipefail
cd "$(dirname "$0")/../man/trace_envelopes"
tmp=$(mktemp -d); trap 'rm -rf "$tmp"' EXIT
mkdir -p "$tmp/figures" "$tmp/generated"
cp trace_inequalities_contrast.tex "$tmp/"
cp ../article/figures/bounds2d_exp1_k2_C1.pdf figures/gardner_contrasts.pdf figures/mc_*.pdf \
   "$tmp/figures/"
cp ../article/generated/trace_characterizations.tex generated/trace_contrast.tex generated/trace_montecarlo.tex \
   generated/trace_montecarlo_figures.tex "$tmp/generated/"
(cd "$tmp" && for i in 1 2 3; do pdflatex -interaction=batchmode -file-line-error trace_inequalities_contrast.tex >/dev/null || true; done)
if grep -q '^\./[A-Za-z_/]*\.tex:[0-9]*:' "$tmp/trace_inequalities_contrast.log" || grep -q 'There were undefined references' "$tmp/trace_inequalities_contrast.log"; then
    grep -E '^\./[A-Za-z_/]*\.tex:[0-9]*:|undefined' "$tmp/trace_inequalities_contrast.log" >&2; echo "arXiv build failed" >&2; exit 1
fi
cp "$tmp/trace_inequalities_contrast.pdf" trace_inequalities_contrast.pdf
rm -f "$tmp"/trace_inequalities_contrast.{aux,log,out,pdf}
tar -czf arxiv_trace_inequalities_contrast.tar.gz -C "$tmp" .
rm -f arxiv_trace_inequalities_contrast.zip && (cd "$tmp" && zip -qr - .) > arxiv_trace_inequalities_contrast.zip
echo "man/trace_envelopes/arxiv_trace_inequalities_contrast.zip: $(tar -tzf arxiv_trace_inequalities_contrast.tar.gz | grep -vc '/$') files; trace_inequalities_contrast.pdf: $(pdfinfo trace_inequalities_contrast.pdf | awk '/Pages/{print $2}') pages"
