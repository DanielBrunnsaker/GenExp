import re

# Matches a single atom once we've isolated it
ATOM_PARSER = re.compile(
    r'^\s*([A-Za-z_]\w*)\s*\(\s*(.*)\s*\)\s*$'
)


def parse_clause(clause):
    """
    Parse a Prolog-like clause into:
        head  – string or None
        atoms – list of (name, args_string)
    Handles:
        - nested parentheses inside arguments
        - quoted strings with commas or parentheses
        - optional trailing period
    """
    text = clause.strip()
    if text.endswith('.'):
        text = text[:-1].rstrip()

    # split head/body on top-level ':-'
    head, body = None, text
    depth = 0
    in_quote = None
    escape = False
    for i in range(len(text) - 1):
        ch = text[i]
        if escape:
            escape = False
            continue
        if ch == '\\':
            escape = True
            continue
        if in_quote:
            if ch == in_quote:
                in_quote = None
            continue
        if ch in ("'", '"'):
            in_quote = ch
            continue
        if ch == '(':
            depth += 1
            continue
        if ch == ')':
            depth = max(depth - 1, 0)
            continue
        if ch == ':' and text[i + 1] == '-' and depth == 0:
            head = text[:i].strip()
            body = text[i + 2:].strip()
            break

    # split body into atoms at top-level commas
    atoms = []
    buf = []
    depth = 0
    in_quote = None
    escape = False
    for ch in body:
        if escape:
            buf.append(ch)
            escape = False
            continue
        if ch == '\\':
            buf.append(ch)
            escape = True
            continue
        if in_quote:
            buf.append(ch)
            if ch == in_quote:
                in_quote = None
            continue
        if ch in ("'", '"'):
            buf.append(ch)
            in_quote = ch
            continue
        if ch == '(':
            depth += 1
            buf.append(ch)
            continue
        if ch == ')':
            depth = max(depth - 1, 0)
            buf.append(ch)
            continue
        if ch == ',' and depth == 0:
            atom = ''.join(buf).strip()
            if atom:
                atoms.append(atom)
            buf = []
            continue
        buf.append(ch)
    final = ''.join(buf).strip()
    if final:
        atoms.append(final)

    # parse each atom into (name, args)
    parsed_atoms = []
    for a in atoms:
        m = ATOM_PARSER.match(a)
        if not m:
            raise ValueError(f"Invalid atom syntax: {a}")
        parsed_atoms.append(
            (m.group(1), [a.replace("'", "") for a in m.group(2).split(",")]))
    return head, parsed_atoms
