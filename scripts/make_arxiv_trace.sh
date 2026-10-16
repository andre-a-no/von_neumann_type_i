#!/usr/bin/env bash
# Build the arXiv upload of the paper on trace characterisations, man/trace_envelopes/arxiv.zip (and arxiv.tar.gz):
# main.tex, the figure and the generated numbers, both written by scripts/trace_characterizations.py into man/article.
# The archive is compiled twice in a clean directory first; the script fails if pdflatex reports an error or an
# undefined reference. It also writes the compiled man/trace_envelopes/main.pdf.
set -euo pipefail
cd "$(dirname "$0")/../man/trace_envelopes"
tmp=$(mktemp -d); trap 'rm -rf "$tmp"' EXIT
mkdir -p "$tmp/figures" "$tmp/generated"
cp main.tex "$tmp/"
cp ../article/figures/trace_characterizations.pdf "$tmp/figures/"
cp ../article/generated/trace_characterizations.tex "$tmp/generated/"
(cd "$tmp" && for i in 1 2 3; do pdflatex -interaction=batchmode -file-line-error main.tex >/dev/null || true; done)
if grep -q '^\./[A-Za-z_/]*\.tex:[0-9]*:' "$tmp/main.log" || grep -q 'There were undefined references' "$tmp/main.log"; then
    grep -E '^\./[A-Za-z_/]*\.tex:[0-9]*:|undefined' "$tmp/main.log" >&2; echo "arXiv build failed" >&2; exit 1
fi
cp "$tmp/main.pdf" main.pdf
rm -f "$tmp"/main.{aux,log,out,pdf}
tar -czf arxiv.tar.gz -C "$tmp" .
rm -f arxiv.zip && (cd "$tmp" && zip -qr - .) > arxiv.zip
echo "man/trace_envelopes/arxiv.zip: $(tar -tzf arxiv.tar.gz | grep -vc '/$') files; main.pdf: $(pdfinfo main.pdf | awk '/Pages/{print $2}') pages"
