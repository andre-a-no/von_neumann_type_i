#!/usr/bin/env bash
# Build the arXiv upload man/article/arxiv_vn_algebra_pytorch.tar.gz (and arxiv_vn_algebra_pytorch.zip): vn_algebra_pytorch.tex, the use-case sections, figures and the generated numbers
# (generated/*.tex). The archive is compiled twice in a clean directory first; the script fails if
# pdflatex reports an error or an undefined reference.
set -euo pipefail
cd "$(dirname "$0")/../man/article"
tmp=$(mktemp -d); trap 'rm -rf "$tmp"' EXIT
mkdir -p "$tmp/figures" "$tmp/generated" "$tmp/usecases"
cp vn_algebra_pytorch.tex sec_control.tex "$tmp/"
cp usecases/*.tex "$tmp/usecases/"
cp figures/*.pdf "$tmp/figures/"
cp generated/*.tex "$tmp/generated/"
(cd "$tmp" && for i in 1 2; do pdflatex -interaction=batchmode -file-line-error vn_algebra_pytorch.tex >/dev/null || true; done)
if grep -q '^\./[A-Za-z_/]*\.tex:[0-9]*:' "$tmp/vn_algebra_pytorch.log" || grep -q 'There were undefined references' "$tmp/vn_algebra_pytorch.log"; then
    grep -E '^\./[A-Za-z_/]*\.tex:[0-9]*:|undefined' "$tmp/vn_algebra_pytorch.log" >&2; echo "arXiv build failed" >&2; exit 1
fi
rm -f "$tmp"/vn_algebra_pytorch.{aux,log,out,pdf}
tar -czf arxiv_vn_algebra_pytorch.tar.gz -C "$tmp" .
rm -f arxiv_vn_algebra_pytorch.zip && (cd "$tmp" && zip -qr - .) > arxiv_vn_algebra_pytorch.zip
echo "man/article/arxiv_vn_algebra_pytorch.tar.gz: $(tar -tzf arxiv_vn_algebra_pytorch.tar.gz | grep -vc '/$') files"
