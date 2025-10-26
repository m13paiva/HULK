# **HULK** Changelog
---

## [1.1.0] - 2025-10-26
### Added
- Switched CLI from **argparse** to **Click** with a custom help renderer:
- **Subcommands** to configure defaults saved locally:
  - `hulk trim` with options:
    - `-ws, --window-size <int>`
    - `-mq, --mean-quality <int>`
  - `hulk tximport` with options:
    - `-m, --mode {raw_counts,length_scaled_tpm,scaled_tpm,dtu_scaled_tpm}`
    - `--ignore-tx-version`
  - Subcommand settings are persisted in a JSON file (default `./.hulk.json`, overridable via `HULK_CONFIG`) and automatically applied when running `hulk ...`.

---

## [1.0.0] - 2025-10-21
### Added
- Initial public release of **HULK**:
---

[1.1.0]: https://github.com/m13paiva/hulk/releases/tag/v1.1.0
[1.0.0]: https://github.com/m13paiva/hulk/releases/tag/v1.0.0

