#!/usr/bin/env bash
set -euo pipefail

# ---- configurable defaults (can be overridden when running installer) ----
IMAGE="${HULK_IMAGE=m13paiva/hulk:latest}"   # <- set HULK_IMAGE=... to customize
APP_NAME="hulk"

# ---- detect container engine ----
ENGINE=""
if command -v docker >/dev/null 2>&1; then
  ENGINE="docker"
elif command -v podman >/dev/null 2>&1; then
  ENGINE="podman"
else
  echo "❌ Neither Docker nor Podman found. Please install one and re-run." >&2
  exit 1
fi

# ---- decide install dir for wrapper ----
TARGET_DIR="/usr/local/bin"
if [ ! -w "$TARGET_DIR" ]; then
  TARGET_DIR="${HOME}/.local/bin"
  mkdir -p "$TARGET_DIR"
  case ":$PATH:" in
    *":${TARGET_DIR}:"*) ;;
    *) echo "ℹ️ Adding ${TARGET_DIR} to PATH in your shell profile."
       SHELLRC="${HOME}/.bashrc"
       [ -n "${ZSH_VERSION:-}" ] && SHELLRC="${HOME}/.zshrc"
       if ! grep -qs "export PATH=.*${TARGET_DIR}" "$SHELLRC"; then
         echo "export PATH=\"${TARGET_DIR}:\$PATH\"" >> "$SHELLRC"
         echo "→ Restart your shell or 'source ${SHELLRC}' after install."
       fi
       ;;
  esac
fi

# ---- local config (persisted on host) ----
CONF_DIR="${HOME}/.hulk"
mkdir -p "${CONF_DIR}/cache" "${CONF_DIR}/ncbi" "${CONF_DIR}"

# Write a sane MultiQC config if missing
MQC_CFG="${CONF_DIR}/multiqc_config.yaml"
if [ ! -f "$MQC_CFG" ]; then
  cat >"$MQC_CFG" <<'YAML'
# multiqc_config.yaml
use_filename_as_sample_name:
  - kallisto
  - fastp
extra_fn_clean_exts:
  - type: regex_keep
    pattern: "(SRR\\d+)"
    module:
      - kallisto
      - fastp
YAML
fi

# ---- pre-pull the image (optional) ----
echo "👉 Pulling image: ${IMAGE} (engine: ${ENGINE})"
${ENGINE} pull "${IMAGE}"

# ---- write wrapper (with installer-resolved image baked in) ----
WRAPPER="${TARGET_DIR}/${APP_NAME}"
cat >"$WRAPPER" <<'WRAP'
#!/usr/bin/env bash
set -euo pipefail

ENGINE="${HULK_ENGINE:-docker}"

# This default gets replaced at install time to match the installer's HULK_IMAGE
IMAGE_DEFAULT="__HULK_IMAGE_DEFAULT__"

# Allow runtime override via env var, otherwise use baked default
IMAGE="${HULK_IMAGE:-$IMAGE_DEFAULT}"

# Resolve engine if overridden / missing
if ! command -v "$ENGINE" >/dev/null 2>&1; then
  if command -v docker >/dev/null 2>&1; then ENGINE="docker"
  elif command -v podman >/dev/null 2>&1; then ENGINE="podman"
  else echo "No container engine found."; exit 1; fi
fi

# Host-side config/cache
CONF_DIR="${HULK_HOME:-$HOME/.hulk}"
mkdir -p "$CONF_DIR" "$CONF_DIR/cache" "$CONF_DIR/ncbi"

# TTY only if interactive
TTY_FLAGS=""
if [ -t 1 ]; then TTY_FLAGS="-it"; fi

# Run as current user so files are owned by you
UIDGID="$(id -u):$(id -g)"

exec "$ENGINE" run --rm $TTY_FLAGS \
  -u "$UIDGID" \
  -e HOME=/work \
  -e MULTIQC_CONFIG_PATH=/config/multiqc_config.yaml \
  -v "$PWD":/app \
  -v "$CONF_DIR":/config \
  -v "$CONF_DIR/cache":/work/.cache \
  -v "$CONF_DIR/ncbi":/work/ncbi \
  -w /app \
  "$IMAGE" "$@"
WRAP
chmod +x "$WRAPPER"

# Replace the placeholder with the resolved image (portable sed)
ESC_IMAGE=$(printf '%s' "$IMAGE" | sed -e 's/[\/&]/\\&/g')
# BSD/GNU sed both accept -i'' (no backup file)
sed -i'' -e "s/__HULK_IMAGE_DEFAULT__/${ESC_IMAGE}/g" "$WRAPPER"

echo "✅ Installed wrapper: $WRAPPER"
echo "   Baked-in image:    $IMAGE"
echo "   Config dir:        $CONF_DIR"
echo
echo "Try:  hulk -h"
