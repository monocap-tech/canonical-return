#!/usr/bin/env python3
"""Execute exactly one staged F audit in a previously assembled Lean checkout.
Usage: python audit_f_stage.py F6_STAGED_AUDIT.json STAGE EVIDENCE_DIR
Final identity requires separately admitted E; operator must enforce that gate.
"""
import hashlib, json, pathlib, re, subprocess, sys

PLAN_BLOB = "866d40011403d87f1f555b19906f31601d3c528e"

def digest(b):
    return hashlib.sha1(b"blob " + str(len(b)).encode() + b"\0" + b).hexdigest()

def main():
    assert len(sys.argv) == 4, __doc__
    path, stage, output = sys.argv[1:]
    raw = pathlib.Path(path).read_bytes()
    assert digest(raw) == PLAN_BLOB, "audit plan drift"
    plan = json.loads(raw)
    selected = [s for s in plan["stages"] if s["stage"] == stage]
    assert len(selected) == 1
    selected = selected[0]
    output = pathlib.Path(output)
    output.mkdir(parents=True, exist_ok=False)
    allowed = set(plan["allowed_axioms"])
    assert allowed == {"propext", "Classical.choice", "Quot.sound"}
    closure = json.loads(pathlib.Path("return109f-closure.json").read_text())
    assert closure["named_targets"] == 111
    results = []
    try:
        for entry in selected["modules"]:
            module = entry["module"]
            source = pathlib.Path(module.replace(".", "/") + ".lean")
            assert digest(source.read_bytes()) == entry["blob"], module
            assert closure["staged"][str(source)]["blob"] == entry["blob"], module
            names = re.findall(r"^#print axioms (.+)$", source.read_text(), re.M)
            assert names == entry["targets"], module
            for suffix, cmd in [
                ("build", ["lake", "build", module]),
                ("audit", ["lake", "env", "lean", str(source)])]:
                with (output / (module + "." + suffix + ".log")).open("wb") as log:
                    p = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT)
                assert p.returncode == 0, (module, suffix, p.returncode)
            log = (output / (module + ".audit.log")).read_text()
            for name in names:
                matches = re.findall(r"'" + re.escape(name) + r"' depends on axioms: \[([^\]]*)\]", log)
                noax = log.count(f"'{name}' does not depend on any axioms")
                assert len(matches) + noax == 1, (name, "missing or duplicate audit")
                axioms = {x.strip() for x in matches[0].split(",") if x.strip()} if matches else set()
                assert axioms <= allowed, (name, sorted(axioms-allowed))
                results.append({"name": name, "axioms": sorted(axioms)})
            print("MODULE PASS", module, len(names), flush=True)
        assert len(results) == selected["count"]
        print("STAGE PASS", stage, len(results), flush=True)
    finally:
        (output / "receipt.json").write_text(json.dumps({
            "stage": stage, "expected": selected["count"], "audited": results,
            "complete": len(results) == selected["count"]
        }, indent=2, sort_keys=True)+"\n")

if __name__ == "__main__":
    main()
