#!/usr/bin/env bash
# Build the arXiv upload of the paper on the boundary time crystal under local dissipation and symmetry-resolved
# entanglement in the XXZ chain, man/btc_sre/arxiv.zip (and arxiv.tar.gz): main.tex, the generated numbers and the
# figures (from scripts/permutational.py and scripts/symmetry_resolved.py in man/article).
# The archive is compiled three times in a clean directory first; the script fails if pdflatex reports an error, an
# undefined reference or citation, or an overfull box wider than 5pt. It also writes the compiled man/btc_sre/main.pdf.
set -euo pipefail
cd "$(dirname "$0")/../man/btc_sre"
tmp=$(mktemp -d); trap 'rm -rf "$tmp"' EXIT
mkdir -p "$tmp/figures" "$tmp/generated"
cp main.tex "$tmp/"
cp ../article/figures/permutational.pdf ../article/figures/symmetry_resolved.pdf "$tmp/figures/"
cp ../article/generated/permutational.tex ../article/generated/symmetry_resolved.tex "$tmp/generated/"
(cd "$tmp" && for i in 1 2 3; do pdflatex -interaction=batchmode -file-line-error main.tex >/dev/null || true; done)
fail=0
if [ ! -f "$tmp/main.pdf" ]; then echo "no main.pdf produced" >&2; fail=1; fi
if grep -q '^\./[A-Za-z_/]*\.tex:[0-9]*:' "$tmp/main.log" || grep -qE 'There were undefined references|Citation .* undefined|Reference .* undefined' "$tmp/main.log"; then
    grep -E '^\./[A-Za-z_/]*\.tex:[0-9]*:|undefined' "$tmp/main.log" >&2; fail=1
fi
wide=$(grep -oE 'Overfull \\[hv]box \(([0-9.]+)pt too' "$tmp/main.log" | grep -oE '[0-9.]+' | awk '$1>5' || true)
if [ -n "$wide" ]; then grep -E 'Overfull' "$tmp/main.log" >&2; fail=1; fi
if [ "$fail" -ne 0 ]; then echo "arXiv build failed" >&2; exit 1; fi
cp "$tmp/main.pdf" main.pdf
rm -f "$tmp"/main.{aux,log,out,pdf}
tar -czf arxiv.tar.gz -C "$tmp" .
rm -f arxiv.zip && (cd "$tmp" && zip -qr - .) > arxiv.zip
echo "man/btc_sre/arxiv.zip: $(tar -tzf arxiv.tar.gz | grep -vc '/$') files; main.pdf: $(pdfinfo main.pdf | awk '/Pages/{print $2}') pages"
