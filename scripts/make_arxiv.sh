#!/usr/bin/env bash
# Build the arXiv upload man/article/arxiv_vn_algebra_pytorch.tar.gz (and arxiv_vn_algebra_pytorch.zip): one
# self-contained vn_algebra_pytorch.tex (sections, use cases and the generated numbers inlined by
# scripts/flatten_arxiv.py, fallbacks resolved) and the figures it uses, all in one directory without
# subdirectories. The archive is compiled twice in a clean directory first; the script fails if pdflatex reports
# an error, an undefined reference or a missing figure.
set -euo pipefail
root="$(cd "$(dirname "$0")/.." && pwd)"
cd "$root/man/article"
tmp=$(mktemp -d); trap 'rm -rf "$tmp"' EXIT
python "$root/scripts/flatten_arxiv.py" . "$tmp/vn_algebra_pytorch.tex"
for f in $(grep -o '\\fig{[^}#]*}' "$tmp/vn_algebra_pytorch.tex" | sed 's/\\fig{\(.*\)}/\1/' | sort -u); do
    cp "figures/$f" "$tmp/"
done
(cd "$tmp" && for i in 1 2; do pdflatex -interaction=batchmode -file-line-error vn_algebra_pytorch.tex >/dev/null || true; done)
if grep -q '^\./[A-Za-z_/]*\.tex:[0-9]*:' "$tmp/vn_algebra_pytorch.log" || grep -q 'not found' "$tmp/vn_algebra_pytorch.log" || grep -q 'There were undefined references' "$tmp/vn_algebra_pytorch.log" \
   || ! grep -q 'Output written' "$tmp/vn_algebra_pytorch.log"; then
    grep -E '^\./[A-Za-z_/]*\.tex:[0-9]*:|undefined' "$tmp/vn_algebra_pytorch.log" >&2; echo "arXiv build failed" >&2; exit 1
fi
echo "arXiv build: $(grep -o '([0-9]* pages' "$tmp/vn_algebra_pytorch.log" | tr -d '(')"
rm -f "$tmp"/vn_algebra_pytorch.{aux,log,out,pdf}
(cd "$tmp" && tar -czf - *) > arxiv_vn_algebra_pytorch.tar.gz
rm -f arxiv_vn_algebra_pytorch.zip && (cd "$tmp" && zip -q - *) > arxiv_vn_algebra_pytorch.zip
echo "man/article/arxiv_vn_algebra_pytorch.tar.gz: $(tar -tzf arxiv_vn_algebra_pytorch.tar.gz | grep -vc '/$') files"
