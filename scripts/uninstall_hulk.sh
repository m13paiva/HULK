#!/usr/bin/env bash
set -euo pipefail

# Confirm
if [[ "${HULK_NUKE:-}" != "1" ]]; then
  echo "This will remove ALL containers/images and uninstall the hulk wrapper & config."
  read -r -p "Type 'YES' to continue: " ANS
  [[ "$ANS" == "YES" ]] || { echo "Aborted."; exit 1; }
fi

# Detect container engine
ENGINE=""
if command -v docker >/dev/null 2>&1; then
  ENGINE="docker"
elif command -v podman >/dev/null 2>&1; then
  ENGINE="podman"
else
  echo "No docker or podman found. Skipping container cleanup."
fi

# Stop & remove ALL containers
if [[ -n "$ENGINE" ]]; then
  echo "Stopping all containers..."
  $ENGINE ps -aq | xargs -r $ENGINE stop

  echo "Removing all containers..."
  $ENGINE ps -aq | xargs -r $ENGINE rm -f -v || true

  echo "Removing all images..."
  $ENGINE images -aq | xargs -r $ENGINE rmi -f || true

  echo "Pruning dangling data (volumes, networks, cache)..."
  # Safe even if nothing to prune
  $ENGINE volume prune -f || true
  if [[ "$ENGINE" == "docker" ]]; then
    docker builder prune -af || true
  fi
  $ENGINE system prune -af --volumes || true
fi

# Remove wrapper(s)
WRAPPERS=(/usr/local/bin/hulk "$HOME/.local/bin/hulk")
for w in "${WRAPPERS[@]}"; do
  if [[ -f "$w" ]]; then
    echo "Removing wrapper: $w"
    rm -f -- "$w"
  fi
done

# Remove config dir (multiqc config, caches)
CONF_DIR="$HOME/.hulk"
if [[ -d "$CONF_DIR" ]]; then
  echo "Removing config dir: $CONF_DIR"
  rm -rf -- "$CONF_DIR"
fi

# Clean PATH/export lines added by installer (best-effort)
clean_rc () {
  local rc="$1"
  [[ -f "$rc" ]] || return 0
  # Remove the exact PATH line the installer may have added
  sed -i '/export PATH=".*\.local\/bin:.*\$PATH"/d' "$rc" || true
  # Remove optional environment overrides
  sed -i '/export HULK_IMAGE=/d' "$rc" || true
  sed -i '/export HULK_ENGINE=/d' "$rc" || true
}

clean_rc "$HOME/.bashrc"
clean_rc "$HOME/.zshrc"

echo "✅ Uninstall complete."
echo "If you modified other shell profiles, you may need to restart your shell."
