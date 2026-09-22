import numpy as np
import math
import time

PHI = (1+math.sqrt(5))/2
FIB = [1,1,2,3,5,8,13,21,34,55,89,144]

# ============================================================
# 3-SAT solver with bitmask + Fibonacci heuristic
# ============================================================

def gen_3sat(n_vars, n_clauses, rng):
    """Random 3-SAT with m/n around phase transition."""
    clauses = []
    for _ in range(n_clauses):
        vars_chosen = rng.choice(n_vars, size=3, replace=False)
        signs = rng.choice([1, -1], size=3)
        clause = [(int(v), int(s)) for v, s in zip(vars_chosen, signs)]
        clauses.append(clause)
    return clauses

def eval_clause(clause, assign):
    """True if at least one literal satisfied."""
    for v, s in clause:
        if v not in assign:
            continue
        if (assign[v] == 1 and s == 1) or (assign[v] == 0 and s == -1):
            return True
    return False

def dpll(clauses, assign, n_vars, max_steps=1_000_000, order_fn=None, counter=None):
    """DPLL with optional variable ordering heuristic."""
    if counter is None:
        counter = [0]
    counter[0] += 1
    if counter[0] > max_steps:
        return None
    
    # Simplify
    new_clauses = []
    for clause in clauses:
        satisfied = False
        rest = []
        for v, s in clause:
            if v in assign:
                if (assign[v] == 1 and s == 1) or (assign[v] == 0 and s == -1):
                    satisfied = True
                    break
            else:
                rest.append((v, s))
        if satisfied:
            continue
        if not rest:
            return False
        new_clauses.append(rest)
    
    if not new_clauses:
        return True
    
    # Unit propagation
    for clause in new_clauses:
        if len(clause) == 1:
            v, s = clause[0]
            val = 1 if s == 1 else 0
            new_assign = dict(assign)
            new_assign[v] = val
            return dpll(new_clauses, new_assign, n_vars, max_steps, order_fn, counter)
    
    # Pure literal elimination
    signs_seen = {}
    for clause in new_clauses:
        for v, s in clause:
            if v in signs_seen:
                if signs_seen[v] != s:
                    signs_seen[v] = 0
            else:
                signs_seen[v] = s
    for v, s in signs_seen.items():
        if s != 0:
            new_assign = dict(assign)
            new_assign[v] = 1 if s == 1 else 0
            return dpll(new_clauses, new_assign, n_vars, max_steps, order_fn, counter)
    
    # Branch
    if order_fn is not None:
        v = order_fn(new_clauses, n_vars, assign)
    else:
        v = new_clauses[0][0][0]
    
    for val in (1, 0):
        new_assign = dict(assign)
        new_assign[v] = val
        result = dpll(new_clauses, new_assign, n_vars, max_steps, order_fn, counter)
        if result is True:
            return True
        if result is None:
            return None
    return False

# ============================================================
# Variable orderings
# ============================================================

def order_fibonacci(clauses, n_vars, assign):
    """Choose variable by Fibonacci index (unassigned)."""
    idx = 0
    while idx < len(FIB):
        v = (FIB[idx] - 1) % n_vars
        if v not in assign:
            return v
        idx += 1
    # Fallback
    for v in range(n_vars):
        if v not in assign:
            return v
    return 0

def order_first(clauses, n_vars, assign):
    """First unassigned variable (default DPLL)."""
    for v in range(n_vars):
        if v not in assign:
            return v
    return 0

def order_most_frequent(clauses, n_vars, assign):
    """Most frequent variable in remaining clauses."""
    counts = {}
    for clause in clauses:
        for v, _ in clause:
            if v not in assign:
                counts[v] = counts.get(v, 0) + 1
    if not counts:
        return 0
    return max(counts, key=counts.get)

# ============================================================
# Main test
# ============================================================
print("="*72)
print("P vs NP: Fibonacci-guided SAT vs standard DPLL")
print("="*72)
print()
print("Ratio m/n = 4.267 (phase transition)")
print()

rng = np.random.default_rng(42)

n_vars_list = [10, 15, 20, 25]
ratio = 4.267
n_formulas = 10
max_steps = 500_000

print(f"{'n':>4} {'m':>5} {'SAT%':>6} | {'first':>12} {'freq':>12} {'fib':>12} | {'fib/first':>10} {'fib/freq':>10}")
print("-"*88)

results = []
for n in n_vars_list:
    m = int(ratio * n)
    n_sat = 0
    steps_first = []
    steps_freq = []
    steps_fib = []
    
    for trial in range(n_formulas):
        clauses = gen_3sat(n, m, rng)
        
        # DPLL with first unassigned
        counter = [0]
        r1 = dpll(clauses, {}, n, max_steps, order_first, counter)
        if r1 is None:
            steps_first.append(max_steps)
        else:
            steps_first.append(counter[0])
            if r1:
                n_sat += 1
        
        # DPLL with most-frequent
        counter = [0]
        r2 = dpll(clauses, {}, n, max_steps, order_most_frequent, counter)
        steps_freq.append(counter[0] if r2 is not None else max_steps)
        
        # DPLL with Fibonacci
        counter = [0]
        r3 = dpll(clauses, {}, n, max_steps, order_fibonacci, counter)
        steps_fib.append(counter[0] if r3 is not None else max_steps)
    
    sf = np.mean(steps_first)
    sq = np.mean(steps_freq)
    sfib = np.mean(steps_fib)
    
    ratio_fib_first = sfib / sf if sf > 0 else 0
    ratio_fib_freq = sfib / sq if sq > 0 else 0
    
    print(f"{n:>4} {m:>5} {n_sat*100//n_formulas:>5}% | {sf:>12.1f} {sq:>12.1f} {sfib:>12.1f} | {ratio_fib_first:>10.3f} {ratio_fib_freq:>10.3f}")
    
    results.append({
        'n': n, 'm': m,
        'sat_rate': n_sat / n_formulas,
        'steps_first': sf, 'steps_freq': sq, 'steps_fib': sfib,
        'ratio_fib_first': ratio_fib_first,
        'ratio_fib_freq': ratio_fib_freq
    })

print()
print("="*72)
print("INTERPRETATION")
print("="*72)
print()
print("ratio_fib_first < 1.0  -> Fibonacci ordering helps")
print("ratio_fib_first ~ 1.0  -> no effect (comparable to baseline)")
print("ratio_fib_first > 1.0  -> Fibonacci ordering hurts")
print()
print("Expected: ratio ~ 1.0 (no exponential speedup).")
print("Fibonacci sequence does NOT give polynomial SAT algorithm.")
print()

# Summary
avg_ratio = np.mean([r['ratio_fib_first'] for r in results])
print(f"Average ratio fib/first = {avg_ratio:.4f}")
if avg_ratio < 0.7:
    print("-> Fibonacci ordering shows some speedup")
elif avg_ratio < 1.3:
    print("-> Fibonacci ordering is neutral (no benefit)")
else:
    print("-> Fibonacci ordering hurts (worse than baseline)")

print()
print("="*72)
print("THEORETICAL NOTE")
print("="*72)
print()
print("P vs NP is a theorem-level problem in complexity theory.")
print("No heuristic ordering of DPLL can resolve it.")
print("This test confirms: Fibonacci sequence is NOT a magic solver.")
print("SAT remains in NP-complete class.")
