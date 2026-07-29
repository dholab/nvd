#!/bin/sh
set -eu

version=3.4.1
private="$PREFIX/libexec/nvd-sra-tools"

rm -rf "$private"
rm -f "$PREFIX/bin/prefetch" "$PREFIX/bin/vdb-validate" "$PREFIX/bin/fastq-dump"
mkdir -p "$private/bin" "$private/lib" "$PREFIX/bin"

case "$target_platform" in
  osx-arm64|osx-64)
    cp "$SRC_DIR/bin/prefetch-orig.$version" "$private/bin/prefetch-orig.$version"
    cp "$SRC_DIR/bin/fastq-dump-orig.$version" "$private/bin/fastq-dump-orig.$version"
    cp "$SRC_DIR/bin/vdb-validate.$version" "$private/bin/vdb-validate.$version"
    cp -R "$SRC_DIR/bin/ncbi" "$private/bin/ncbi"
    cp -R "$SRC_DIR/schema" "$private/schema"
    ;;

  linux-64|linux-aarch64)
    case "$target_platform" in
      linux-64)
        manifest=sha256:1bf9aa259adc315df1cf9d6f8e9d5de7146f1d167e4b7abdc7a61438423da795
        loader=ld-musl-x86_64.so.1
        ;;
      linux-aarch64)
        manifest=sha256:5ad41a0068a1b9b1c3fcba54ec6cc8605760a379b231df82167c4ea70f166760
        loader=ld-musl-aarch64.so.1
        ;;
    esac

    image="$BUILD_DIR/sra-tools-image"
    rootfs="$BUILD_DIR/sra-tools-rootfs"
    rm -rf "$image" "$rootfs"
    mkdir -p "$image" "$rootfs"

    skopeo --registries-conf /dev/null copy \
      "docker://docker.io/ncbi/sra-tools@$manifest" "dir:$image"

    jq -r '.layers[].digest' "$image/manifest.json" |
      while IFS= read -r digest; do
        tar -xf "$image/${digest#sha256:}" -C "$rootfs"
      done

    cp "$rootfs/usr/local/bin/prefetch-orig.$version" "$private/bin/prefetch-orig.$version"
    cp "$rootfs/usr/local/bin/fastq-dump-orig.$version" "$private/bin/fastq-dump-orig.$version"
    cp "$rootfs/usr/local/bin/vdb-validate.$version" "$private/bin/vdb-validate.$version"
    cp -R "$rootfs/usr/local/bin/ncbi" "$private/bin/ncbi"
    cp "$rootfs/lib/$loader" "$private/lib/$loader"
    ln -s "$loader" "$private/lib/ld-musl.so.1"
    ;;

  *)
    printf 'unsupported target platform: %s\n' "$target_platform" >&2
    exit 1
    ;;
esac

chmod 755 \
  "$private/bin/prefetch-orig.$version" \
  "$private/bin/fastq-dump-orig.$version" \
  "$private/bin/vdb-validate.$version"

for tool in prefetch vdb-validate fastq-dump; do
  cp "$RECIPE_DIR/launcher" "$PREFIX/bin/$tool"
  chmod 755 "$PREFIX/bin/$tool"
done
