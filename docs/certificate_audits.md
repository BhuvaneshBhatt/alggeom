# Certificate-field audits

`tests/certificate_fields.json` classifies the stored fields of each registered certificate family. New or renamed fields fail the audit until their role is recorded. No field is skipped merely because it is named method, source or modulus.

Proof fields have explicit mutations and replay must reject them. Non-proof fields require a reason. semialg additionally records conditional branches as inactive when absent, constructor/schema rejections separately from exact replay, and per-root overrides when a shared component type plays a different role in another certificate. Unexpected verifier crashes fail the tests.

Set `CERTIFICATE_AUDIT_REPORT=certificate-audit.jsonl` during pytest to retain field paths and outcomes. The report distinguishes mutated evidence, provenance, diagnostics and inactive branches. This is semantic certificate mutation coverage, not a source-code mutation score. The registry covers the exercised families; it does not establish that every possible malformed Python object or every certificate family in the package has been enumerated.

Independent arithmetic identities, coordinate transports and prescribed fibers complement replay tests so a constructor and verifier cannot serve as each other's only oracle. See the package testing guide for the registered workflows.
