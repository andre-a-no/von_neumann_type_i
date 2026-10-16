#!/usr/bin/env bash
# Build the arXiv upload of the paper on the Connes--Stormer entropy of two MASAs, man/connes_stormer/arxiv_connes_stormer_masas.zip (and
# arxiv_connes_stormer_masas.tar.gz): connes_stormer_masas.tex, the generated numbers and the figure (both from scripts/connes_stormer.py in man/article).
# The archive is compiled three times in a clean directory first; the script fails if pdflatex reports an error or an
# undefined reference. It also writes the compiled man/connes_stormer/connes_stormer_masas.pdf.
set -euo pipefail
cd "$(dirname "$0")/../man/connes_stormer"
tmp=$(mktemp -d); trap 'rm -rf "$tmp"' EXIT
mkdir -p "$tmp/figures" "$tmp/generated"
cp connes_stormer_masas.tex "$tmp/"
cp ../article/figures/connes_stormer.pdf "$tmp/figures/"
cp ../article/generated/connes_stormer.tex "$tmp/generated/"
(cd "$tmp" && for i in 1 2 3; do pdflatex -interaction=batchmode -file-line-error connes_stormer_masas.tex >/dev/null || true; done)
if [ ! -f "$tmp/connes_stormer_masas.pdf" ] || grep -q '^\./[A-Za-z_/]*\.tex:[0-9]*:' "$tmp/connes_stormer_masas.log" \
        || grep -q -E 'There were undefined references|Reference .* undefined|Citation .* undefined' "$tmp/connes_stormer_masas.log"; then
    grep -E '^\./[A-Za-z_/]*\.tex:[0-9]*:|undefined' "$tmp/connes_stormer_masas.log" >&2 || true; echo "arXiv build failed" >&2; exit 1
fi
awk '/^Overfull \\hbox \(/ { w = $3; sub(/^\(/, "", w); sub(/pt$/, "", w); if (w + 0 > 5) print "warning: " $0 }' "$tmp/connes_stormer_masas.log" >&2
cp "$tmp/connes_stormer_masas.pdf" connes_stormer_masas.pdf
rm -f "$tmp"/connes_stormer_masas.{aux,log,out,pdf}
tar -czf arxiv_connes_stormer_masas.tar.gz -C "$tmp" .
rm -f arxiv_connes_stormer_masas.zip && (cd "$tmp" && zip -qr - .) > arxiv_connes_stormer_masas.zip
echo "man/connes_stormer/arxiv_connes_stormer_masas.zip: $(tar -tzf arxiv_connes_stormer_masas.tar.gz | grep -vc '/$') files; connes_stormer_masas.pdf: $(pdfinfo connes_stormer_masas.pdf | awk '/Pages/{print $2}') pages"
