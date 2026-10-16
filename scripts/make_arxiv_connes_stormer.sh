#!/usr/bin/env bash
# Build the arXiv upload of the paper on the Connes--Stormer entropy of two MASAs, man/connes_stormer/arxiv.zip (and
# arxiv.tar.gz): main.tex, the generated numbers and the figure (both from scripts/connes_stormer.py in man/article).
# The archive is compiled three times in a clean directory first; the script fails if pdflatex reports an error or an
# undefined reference. It also writes the compiled man/connes_stormer/main.pdf.
set -euo pipefail
cd "$(dirname "$0")/../man/connes_stormer"
tmp=$(mktemp -d); trap 'rm -rf "$tmp"' EXIT
mkdir -p "$tmp/figures" "$tmp/generated"
cp main.tex "$tmp/"
cp ../article/figures/connes_stormer.pdf "$tmp/figures/"
cp ../article/generated/connes_stormer.tex "$tmp/generated/"
(cd "$tmp" && for i in 1 2 3; do pdflatex -interaction=batchmode -file-line-error main.tex >/dev/null || true; done)
if [ ! -f "$tmp/main.pdf" ] || grep -q '^\./[A-Za-z_/]*\.tex:[0-9]*:' "$tmp/main.log" \
        || grep -q -E 'There were undefined references|Reference .* undefined|Citation .* undefined' "$tmp/main.log"; then
    grep -E '^\./[A-Za-z_/]*\.tex:[0-9]*:|undefined' "$tmp/main.log" >&2 || true; echo "arXiv build failed" >&2; exit 1
fi
awk '/^Overfull \\hbox \(/ { w = $3; sub(/^\(/, "", w); sub(/pt$/, "", w); if (w + 0 > 5) print "warning: " $0 }' "$tmp/main.log" >&2
cp "$tmp/main.pdf" main.pdf
rm -f "$tmp"/main.{aux,log,out,pdf}
tar -czf arxiv.tar.gz -C "$tmp" .
rm -f arxiv.zip && (cd "$tmp" && zip -qr - .) > arxiv.zip
echo "man/connes_stormer/arxiv.zip: $(tar -tzf arxiv.tar.gz | grep -vc '/$') files; main.pdf: $(pdfinfo main.pdf | awk '/Pages/{print $2}') pages"
