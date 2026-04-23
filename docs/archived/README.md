# BioAssert Project — Export Package

**Package / repo:** `bioassert`
**Dataset release:** `BioAssert-NSCLC-v1`

Three documents. Read in this order:

1. **conversation_context.md** — the reasoning chain. Why this project, what problem it solves, what makes it methodologically defensible. Read first.

2. **project_spec.md** — the full architectural specification. Seven-layer architecture, ontology schema, complexity stratification, evaluation harness, repo layout, development phases, non-negotiables. This is the reference doc to paste into every Claude Code session.

3. **claude_code_kickoff.md** — the exact prompt to start the first development session, plus the validation checklist for Phase 1 output.

## Next action

Open Claude Code in a new repo directory. Attach all three markdown files as context. Paste the prompt from claude_code_kickoff.md.

Keep each Claude Code session tightly scoped to one architectural layer. Don't try to build L1 through L7 in one go.
