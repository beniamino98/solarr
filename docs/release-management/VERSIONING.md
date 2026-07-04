# Release-managed app Release Management

Markdown-only release register for the Release-managed app Shiny app.
This file is intended to be machine-readable by the release-management dashboard.
Closed issues are recorded here only after verification evidence exists.

## Current Version

| Field | Value |
| --- | --- |
| App name | `Release-managed app` |
| App prefix | `SOLARR` |
| Issue prefix | `PKG-SOLARR` |
| Current production version | `0.0.0` |
| Active / latest release branch | `SOLARR-0.0.1` |
| Target release version | `0.0.1` |
| Rules file | `docs/release-management/RULES.md` |
| Tracker | `docs/release-management/TRACKER.md` |
| Versioning rule | `Only closed and verified issues are listed here. Open issues must stay in TRACKER.md.` |

## Release Register

### `SOLARR-0.0.1`

Target production version after merge: `0.0.1`  
Baseline production version: `0.0.0`  
Branch: `SOLARR-0.0.1`  
Merge date: ``  
Release type: `PATCH`

| Commit | Date closed | Change type | Issues | Opened in version | Summary | Verification | Production version after merge |
| --- | --- | --- | --- | --- | --- | --- | --- |

## Field Guide

- `Commit`: commit identifier/message prefix in the form `PREFIX-TARGET_VERSION-###`.
- `Issues`: one or more issue IDs resolved by the commit, separated by semicolons when grouped.
- `Opened in version`: production version where the issue was originally registered.
- `Production version after merge`: production version reached when the branch is merged.

## Editing Notes

Do not rename, remove, or reorder the release register columns unless the Shiny parser is updated at the same time.
Operational versioning, branch, commit, push, and merge rules live in `RULES.md`.
Open issues must stay in `TRACKER.md` until verified and closed.
