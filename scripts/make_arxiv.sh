#!/usr/bin/env bash
# Build the arXiv upload man/article/arxiv.tar.gz: main.tex, the use-case sections, figures, listings and the generated numbers
# (generated/*.tex). The archive is compiled twice in a clean directory first; the script fails if
# pdflatex reports an error or an undefined reference.
set -euo pipefail
cd "$(dirname "$0")/../man/article"
tmp=$(mktemp -d); trap 'rm -rf "$tmp"' EXIT
mkdir -p "$tmp/figures" "$tmp/listings" "$tmp/generated" "$tmp/usecases"
cp main.tex "$tmp/"
cp usecases/*.tex "$tmp/usecases/"
cp figures/*.pdf "$tmp/figures/"
cp listings/*.py "$tmp/listings/"
cp generated/*.tex "$tmp/generated/"
(cd "$tmp" && for i in 1 2; do pdflatex -interaction=batchmode -file-line-error main.tex >/dev/null || true; done)
if grep -q '^\./[A-Za-z_/]*\.tex:[0-9]*:' "$tmp/main.log" || grep -q 'There were undefined references' "$tmp/main.log"; then
    grep -E '^\./[A-Za-z_/]*\.tex:[0-9]*:|undefined' "$tmp/main.log" >&2; echo "arXiv build failed" >&2; exit 1
fi
rm -f "$tmp"/main.{aux,log,out,pdf}
tar -czf arxiv.tar.gz -C "$tmp" .
echo "man/article/arxiv.tar.gz: $(tar -tzf arxiv.tar.gz | grep -vc '/$') files"
