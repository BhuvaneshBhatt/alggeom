"""Field-complete mutation cases with explicit per-certificate policies."""

import json
import os
from dataclasses import fields, replace
from pathlib import Path

POLICIES = json.loads(
    Path(__file__).with_name("certificate_fields.json").read_text(encoding="utf-8")
)


def audit_fields(certificate, mutations, verify):
    """Require an explicit test or documented classification for every stored field."""
    name = type(certificate).__module__ + "." + type(certificate).__qualname__
    assert name in POLICIES, f"unclassified certificate: {name}"
    policy = POLICIES[name]
    actual = {field.name for field in fields(certificate)}
    assert set(policy) == actual, f"unclassified fields in {name}: {actual - set(policy)}"
    proof = {key for key, value in policy.items() if value["role"] == "proof"}
    assert set(mutations) == proof, f"mutation omissions in {name}: {proof - set(mutations)}"
    assert verify(certificate)
    report = []
    for key, entry in policy.items():
        if key not in proof:
            assert entry.get("reason")
            report.append({"path": key, "status": entry["role"], "reason": entry["reason"]})
            continue
        changed = replace(certificate, **{key: mutations[key]})
        assert not verify(changed), f"altered field accepted: {name}.{key}"
        report.append({"path": key, "status": "replay-rejected"})
    record = {
        "certificate": name,
        "fields": report,
        "scope": "stored root fields; nested identities have separate cases",
    }
    target = os.environ.get("CERTIFICATE_AUDIT_REPORT")
    if target:
        with Path(target).open("a", encoding="utf-8") as stream:
            stream.write(json.dumps(record) + "\n")
    return record
