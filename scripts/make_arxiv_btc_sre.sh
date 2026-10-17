#!/usr/bin/env bash
# Build the arXiv upload of the paper on the boundary time crystal under local dissipation and symmetry-resolved
# entanglement in the XXZ chain, man/btc_sre/arxiv_btc_sre_symmetry_sectors.zip (and arxiv_btc_sre_symmetry_sectors.tar.gz): btc_sre_symmetry_sectors.tex, the generated numbers and the
# figures (from scripts/permutational.py and scripts/symmetry_resolved.py in man/article).
# The archive is compiled three times in a clean directory first; the script fails if pdflatex reports an error, an
# undefined reference or citation, or an overfull box wider than 5pt. It also writes the compiled man/btc_sre/btc_sre_symmetry_sectors.pdf.
set -euo pipefail
cd "$(dirname "$0")/../man/btc_sre"
tmp=$(mktemp -d); trap 'rm -rf "$tmp"' EXIT
mkdir -p "$tmp/figures" "$tmp/generated"
cp btc_sre_symmetry_sectors.tex "$tmp/"
cp ../article/figures/permutational.pdf ../article/figures/symmetry_resolved.pdf "$tmp/figures/"
cp ../article/generated/permutational.tex ../article/generated/symmetry_resolved.tex "$tmp/generated/"
(cd "$tmp" && for i in 1 2 3; do pdflatex -interaction=batchmode -file-line-error btc_sre_symmetry_sectors.tex >/dev/null || true; done)
fail=0
if [ ! -f "$tmp/btc_sre_symmetry_sectors.pdf" ]; then echo "no btc_sre_symmetry_sectors.pdf produced" >&2; fail=1; fi
if grep -q '^\./[A-Za-z_/]*\.tex:[0-9]*:' "$tmp/btc_sre_symmetry_sectors.log" || grep -qE 'There were undefined references|Citation .* undefined|Reference .* undefined' "$tmp/btc_sre_symmetry_sectors.log"; then
    grep -E '^\./[A-Za-z_/]*\.tex:[0-9]*:|undefined' "$tmp/btc_sre_symmetry_sectors.log" >&2; fail=1
fi
wide=$(grep -oE 'Overfull \\[hv]box \(([0-9.]+)pt too' "$tmp/btc_sre_symmetry_sectors.log" | grep -oE '[0-9.]+' | awk '$1>5' || true)
if [ -n "$wide" ]; then grep -E 'Overfull' "$tmp/btc_sre_symmetry_sectors.log" >&2; fail=1; fi
if [ "$fail" -ne 0 ]; then echo "arXiv build failed" >&2; exit 1; fi
cp "$tmp/btc_sre_symmetry_sectors.pdf" btc_sre_symmetry_sectors.pdf
rm -f "$tmp"/btc_sre_symmetry_sectors.{aux,log,out,pdf}
tar -czf arxiv_btc_sre_symmetry_sectors.tar.gz -C "$tmp" .
rm -f arxiv_btc_sre_symmetry_sectors.zip && (cd "$tmp" && zip -qr - .) > arxiv_btc_sre_symmetry_sectors.zip
echo "man/btc_sre/arxiv_btc_sre_symmetry_sectors.zip: $(tar -tzf arxiv_btc_sre_symmetry_sectors.tar.gz | grep -vc '/$') files; btc_sre_symmetry_sectors.pdf: $(pdfinfo btc_sre_symmetry_sectors.pdf | awk '/Pages/{print $2}') pages"
