# Open Questions Merge: Master Log

**Created:** 2026-04-12
**Scope:** Tracks the per-batch human review of Cowork's proposal files
in `review/OPEN_QUESTIONS_INDEX_PROPOSAL_*.md`. See
`ClaudeTasks/TASK_OPEN_QUESTIONS_MERGE.md` for the protocol.

Proposal files remain untouched (Cowork's original classifications are
preserved as archaeology). This master log and the per-batch logs below
are the human-review layer on top of them.

---

## Batch progress

| Batch | Entries | Done | Log file |
|-------|---------|------|----------|
| methodology | 8 | 2 | `MERGE_LOG_methodology.md` |
| closed-form | 9 | 0 | (not yet started) |
| interpretation | 18 | 0 | (not yet started) |
| scope-extension | 26 | 0 | (not yet started) |
| hardware-test | 38 | 0 | (not yet started) |
| math-proof | 43 | 0 | (not yet started) |
| numerical-verification | 62 | 0 | (not yet started) |
| untagged 8a | ~41 | 0 | (not yet started) |
| untagged 8b | ~41 | 0 | (not yet started) |
| untagged 8c | ~41 | 0 | (not yet started) |

**Total:** 2 / 327

**Last updated:** 2026-04-12

---

## Status vocabulary used in per-batch logs

Beyond Cowork's six classes (`open`, `resolved`, `partially-resolved`,
`superseded`, `obsolete`, `needs-human`), this human-review layer uses
two additional meta-categories introduced in the merge protocol:

- **`todo`**: Self-assignment, not a research question. The scientific
  content is settled; only a code/pipeline/scaffolding change remains.
  Resolution lives outside science docs (in `ClaudeTasks/` or directly
  in code). Cowork cannot reach this status on its own because the
  resolution location is out of scope for its scan.

- **`by-design-open`**: A caveat explicitly declared as a limit of a
  Tier-5 synthesis frame. The gap is known, owned, and not expected to
  be closed within the current framing. Distinct from `open` in that
  it is not a research debt but a deliberately acknowledged boundary.

When an entry is classified as `todo` and then executed in the same or
a later session, its final status is upgraded to `resolved` with a
rationale pointing to the implementing commit.
