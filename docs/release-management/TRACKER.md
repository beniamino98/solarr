# Release-managed app Issue Tracker

Markdown-only tracker for the Release-managed app Shiny app.
This file is intended to be machine-readable by the release-management dashboard.
Do not create or maintain Excel trackers for this app.

## Tracker Metadata

| Field | Value |
| --- | --- |
| App name | `Release-managed app` |
| App prefix | `SOLARR` |
| Issue prefix | `PKG-SOLARR` |
| Current production version | `0.0.0` |
| Active release branch | `SOLARR-0.0.1` |
| Target release version | `0.0.1` |
| Rules file | `docs/release-management/RULES.md` |
| Versioning file | `docs/release-management/VERSIONING.md` |
| Tracker rule | `Only open, in-progress, or blocked issues are listed here. Closed issues must be moved to VERSIONING.md.` |

## Open Issues

| Issue ID | Opened in version | Change type | Title | Severity | Owner module | Affected files/functions | Evidence | Recommended fix | Status | Opened | Verification plan | Notes |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |

## Field Guide

- `Issue ID`: progressive issue identifier independent from app version.
- `Opened in version`: production version where the issue was noticed and registered.
- `Change type`: expected release impact, one of `PATCH`, `MINOR`, or `MAJOR`.
- `Status`: one of `Open`, `In progress`, or `Blocked`; closed issues must be moved to `VERSIONING.md`.
- `Evidence`, `Recommended fix`, `Verification plan`, and `Notes`: editable long text fields.

## Editing Notes

Do not rename, remove, or reorder the `Open Issues` columns unless the Shiny parser is updated at the same time.
Use backticks around IDs, versions, dates, paths, and function names.
