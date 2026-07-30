{
  description = "Reproducible dev shell for the `NVD` bioinformatic processing pipeline";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixpkgs-unstable";
    flake-utils.url = "github:numtide/flake-utils";
  };

  outputs =
    {
      self,
      nixpkgs,
      flake-utils,
      ...
    }:
    flake-utils.lib.eachDefaultSystem (
      system:
      let
        pkgs = import nixpkgs {
          inherit system;
        };

        pixiVersion = "0.74.0";
        pixiRelease = {
          aarch64-darwin = {
            target = "aarch64-apple-darwin";
            hash = "sha256-t8kqwVMXHSXEanOHJyvaGs5tn530b9NLPJc3MAsyJBU=";
          };
          x86_64-darwin = {
            target = "x86_64-apple-darwin";
            hash = "sha256-II5HVbzfrIqWxTAU/WOe/nUy5qswGwNgq03vK+Max98=";
          };
          aarch64-linux = {
            target = "aarch64-unknown-linux-musl";
            hash = "sha256-hJ3dnaP82/yZpZvi8zlzEX0ZrVOjQIYewrmahqm2F98=";
          };
          x86_64-linux = {
            target = "x86_64-unknown-linux-musl";
            hash = "sha256-BuMYXJdAr5/9NFYQHXvw6uw7KUrxgIkijj8lnDneC2Q=";
          };
        }.${system};
        pixi = pkgs.stdenvNoCC.mkDerivation {
          pname = "pixi";
          version = pixiVersion;
          src = pkgs.fetchurl {
            url = "https://github.com/prefix-dev/pixi/releases/download/v${pixiVersion}/pixi-${pixiRelease.target}";
            hash = pixiRelease.hash;
          };
          dontUnpack = true;
          installPhase = ''
            runHook preInstall
            install -Dm755 "$src" "$out/bin/pixi"
            runHook postInstall
          '';
        };

      in
      {
        devShells.default = pkgs.mkShell {
          name = "NVD";

          buildInputs = [
            pkgs.stdenv
            pkgs.gcc
            pkgs.curl
            pkgs.wget
            pkgs.openjdk
            pkgs.git
            pkgs.cmake
            pkgs.libxml2
            pkgs.libxslt
            pkgs.libffi
            pixi
            pkgs.graphviz
          ];

          shellHook = ''
            echo "🔧 Entering NVD dev shell"
            export PS1="(NVD) $PS1"
            if [ ! -d .pixi/envs/default ]; then
              echo "Pixi env not found. Running install..."
              pixi install --frozen
            fi

            export PATH="$PWD/.pixi/envs/default/bin:$PATH"
          '';
        };
      }
    );
}
