"""
Biomedical Data Structures and Algorithms Project (U25ADG02)

Module 1 (Unit 3): Expression tree for clinical formulas (BMI, BSA, MAP)
Module 2 (Unit 4): DFS for protein interaction pathways
Module 3 (Unit 5): Image-filtering operation cost

Team: 25BM050 Srikrishna S | 25BM051 Srinath R |
      25BM052 Srinithi P | 25BM053 Sruthika V

All data used here is synthetic and for education only.
"""

# =====================================================================
# MODULE 1 - EXPRESSION TREE FOR CLINICAL FORMULAS
# =====================================================================

class ExprNode:
    def __init__(self, value, left=None, right=None):
        self.value = value
        self.left = left
        self.right = right


PRECEDENCE = {'+': 1, '-': 1, '*': 2, '/': 2, '^': 3}
RIGHT_ASSOCIATIVE = {'^'}


def tokenize(expression):
    tokens, i = [], 0

    while i < len(expression):
        ch = expression[i]

        if ch.isspace():
            i += 1

        elif ch.isdigit() or ch == '.':
            j = i
            while j < len(expression) and (
                expression[j].isdigit() or expression[j] == '.'
            ):
                j += 1
            tokens.append(expression[i:j])
            i = j

        elif ch.isalpha() or ch == '_':
            j = i
            while j < len(expression) and (
                expression[j].isalnum() or expression[j] == '_'
            ):
                j += 1
            tokens.append(expression[i:j])
            i = j

        elif ch in PRECEDENCE or ch in '()':
            tokens.append(ch)
            i += 1

        else:
            raise ValueError("Invalid character in formula: " + ch)

    return tokens


def is_operand(token):
    return token not in PRECEDENCE and token not in ('(', ')')


def to_postfix(tokens):
    """Convert infix tokens to postfix using the Shunting-yard algorithm."""
    output, stack = [], []

    for token in tokens:
        if is_operand(token):
            output.append(token)

        elif token == '(':
            stack.append(token)

        elif token == ')':
            while stack and stack[-1] != '(':
                output.append(stack.pop())

            if not stack:
                raise ValueError("Mismatched parentheses")

            stack.pop()

        else:
            while (
                stack
                and stack[-1] != '('
                and (
                    PRECEDENCE[stack[-1]] > PRECEDENCE[token]
                    or (
                        PRECEDENCE[stack[-1]] == PRECEDENCE[token]
                        and token not in RIGHT_ASSOCIATIVE
                    )
                )
            ):
                output.append(stack.pop())

            stack.append(token)

    while stack:
        if stack[-1] == '(':
            raise ValueError("Mismatched parentheses")
        output.append(stack.pop())

    return output


def build_tree(postfix):
    stack = []

    for token in postfix:
        if token in PRECEDENCE:
            if len(stack) < 2:
                raise ValueError("Operator is missing an operand")

            right = stack.pop()
            left = stack.pop()
            stack.append(ExprNode(token, left, right))

        else:
            stack.append(ExprNode(token))

    if len(stack) != 1:
        raise ValueError("Malformed formula")

    return stack[0]


def evaluate(node, values):
    """Evaluate an expression tree using supplied variable values."""
    if node.left is None and node.right is None:
        try:
            return float(node.value)
        except ValueError:
            if node.value not in values:
                raise KeyError("No value supplied for " + node.value)
            return float(values[node.value])

    a = evaluate(node.left, values)
    b = evaluate(node.right, values)

    if node.value == '+':
        return a + b

    if node.value == '-':
        return a - b

    if node.value == '*':
        return a * b

    if node.value == '/':
        if b == 0:
            raise ZeroDivisionError("Division by zero in formula")
        return a / b

    if a < 0 and b != int(b):
        raise ValueError("Negative base with fractional power")

    return a ** b


def infix(node):
    if node.left is None:
        return node.value

    return (
        '('
        + infix(node.left)
        + ' '
        + node.value
        + ' '
        + infix(node.right)
        + ')'
    )


def preorder(node):
    if node is None:
        return []

    return [node.value] + preorder(node.left) + preorder(node.right)


def postorder(node):
    if node is None:
        return []

    return postorder(node.left) + postorder(node.right) + [node.value]


def tree_height(node):
    if node is None:
        return 0

    return 1 + max(tree_height(node.left), tree_height(node.right))


def count_nodes(node):
    if node is None:
        return 0

    return 1 + count_nodes(node.left) + count_nodes(node.right)


def print_tree(node, depth=0):
    """Print the expression tree sideways."""
    if node is not None:
        print_tree(node.right, depth + 1)
        print('    ' * depth + str(node.value))
        print_tree(node.left, depth + 1)


FORMULAS = {
    'BMI': ('wt / ht_m ^ 2', 'kg/m^2'),
    'BSA': ('( ht_cm * wt / 3600 ) ^ 0.5', 'm^2'),
    'MAP': ('dbp + ( sbp - dbp ) / 3', 'mmHg'),
}


def make_tree(expression):
    return build_tree(to_postfix(tokenize(expression)))


# =====================================================================
# MODULE 2 - DFS FOR PROTEIN INTERACTION PATHWAYS
# =====================================================================

def build_graph(proteins, interactions):
    graph = {p: [] for p in proteins}

    for a, b in interactions:
        if a not in graph or b not in graph:
            raise ValueError("Interaction references an undefined protein")

        if b not in graph[a]:
            graph[a].append(b)
            graph[b].append(a)

    return graph


def dfs_recursive(graph, node, visited=None, order=None):
    if visited is None:
        visited, order = set(), []

    visited.add(node)
    order.append(node)

    for neighbour in graph[node]:
        if neighbour not in visited:
            dfs_recursive(graph, neighbour, visited, order)

    return order


def dfs_iterative(graph, start):
    visited, order, stack = set(), [], [start]

    while stack:
        node = stack.pop()

        if node in visited:
            continue

        visited.add(node)
        order.append(node)

        for neighbour in reversed(graph[node]):
            if neighbour not in visited:
                stack.append(neighbour)

    return order


def find_pathway(graph, source, target):
    """Find the first pathway discovered by DFS with backtracking."""
    visited, path = set(), []

    def visit(node):
        visited.add(node)
        path.append(node)

        if node == target:
            return True

        for neighbour in graph[node]:
            if neighbour not in visited and visit(neighbour):
                return True

        path.pop()
        return False

    return list(path) if visit(source) else None


def all_pathways(graph, source, target, max_length=6):
    """Find every simple pathway with at most max_length proteins."""
    pathways = []

    def visit(node, path):
        if node == target:
            pathways.append(list(path))
            return

        if len(path) == max_length:
            return

        for neighbour in graph[node]:
            if neighbour not in path:
                path.append(neighbour)
                visit(neighbour, path)
                path.pop()

    visit(source, [source])
    return pathways


def connected_components(graph):
    seen, groups = set(), []

    for protein in graph:
        if protein not in seen:
            group = dfs_iterative(graph, protein)
            seen.update(group)
            groups.append(group)

    return groups


def has_cycle(graph):
    """Detect a cycle in the undirected interaction network."""
    visited = set()

    def visit(node, parent):
        visited.add(node)

        for neighbour in graph[node]:
            if neighbour not in visited:
                if visit(neighbour, node):
                    return True
            elif neighbour != parent:
                return True

        return False

    return any(
        protein not in visited and visit(protein, None)
        for protein in graph
    )


PROTEINS = [
    'TP53', 'MDM2', 'MDM4', 'ATM', 'CHEK2', 'BRCA1',
    'RAD51', 'EP300', 'EGFR', 'GRB2', 'SOS1', 'KRAS'
]

INTERACTIONS = [
    ('TP53', 'MDM2'),
    ('MDM2', 'MDM4'),
    ('TP53', 'ATM'),
    ('TP53', 'EP300'),
    ('ATM', 'CHEK2'),
    ('CHEK2', 'BRCA1'),
    ('BRCA1', 'RAD51'),
    ('TP53', 'BRCA1'),
    ('EGFR', 'GRB2'),
    ('GRB2', 'SOS1'),
    ('SOS1', 'KRAS'),
]


# =====================================================================
# MODULE 3 - IMAGE-FILTER OPERATION COUNT
# =====================================================================

def filter_operation_count(height, width, k):
    if height <= 0 or width <= 0:
        raise ValueError("Image dimensions must be positive")

    if k <= 0 or k % 2 == 0:
        raise ValueError("Filter size must be a positive odd number")

    if k > height or k > width:
        raise ValueError("Filter window cannot exceed image dimensions")

    out_h = height - k + 1
    out_w = width - k + 1

    visits = 0

    for r in range(out_h):
        for c in range(out_w):
            for i in range(k):
                for j in range(k):
                    visits += 1

    return out_h, out_w, visits


# =====================================================================
# MENU-DRIVEN INTERFACE
# =====================================================================

def read_int(prompt, low=None, high=None):
    while True:
        text = input(prompt).strip()

        try:
            value = int(text)
        except ValueError:
            print("Invalid input: please enter a whole number.")
            continue

        if (
            (low is not None and value < low)
            or (high is not None and value > high)
        ):
            print("Invalid input: choose a value between", low, "and", high)
            continue

        return value


def show_tree_report(expression, values):
    tree = make_tree(expression)

    print("Formula (infix)    :", infix(tree))
    print("Preorder           :", ' '.join(preorder(tree)))
    print("Postorder          :", ' '.join(postorder(tree)))
    print(
        "Nodes / height     :",
        count_nodes(tree),
        "/",
        tree_height(tree)
    )
    print("Tree (sideways):")
    print_tree(tree)

    return evaluate(tree, values)


def run_expression_module():
    patient = {
        'wt': 70,
        'ht_m': 1.75,
        'ht_cm': 175,
        'sbp': 120,
        'dbp': 80
    }

    print("\n--- Module 1: Expression tree for clinical formulas ---")
    print("1. BMI   2. BSA (Mosteller)   3. MAP   4. Custom formula")
    print("0. Back")

    choice = read_int("Choice: ", 0, 4)

    try:
        if choice in (1, 2, 3):
            name = ['BMI', 'BSA', 'MAP'][choice - 1]
            expression, unit = FORMULAS[name]

            print("Synthetic patient values:", patient)

            result = show_tree_report(expression, patient)
            print(name, "=", round(result, 2), unit)

        elif choice == 4:
            expression = input(
                "Enter formula (+ - * / ^, parentheses): "
            )

            values = {}

            for token in tokenize(expression):
                if is_operand(token) and token not in values:
                    try:
                        float(token)
                    except ValueError:
                        values[token] = float(
                            input("Value of " + token + ": ")
                        )

            result = show_tree_report(expression, values)
            print("Result =", round(result, 4))

    except (ValueError, KeyError, ZeroDivisionError) as error:
        print("Error:", error)


def run_protein_module():
    graph = build_graph(PROTEINS, INTERACTIONS)

    print("\n--- Module 2: DFS for protein interaction pathways ---")
    print("1. DFS traversal 2. First pathway 3. All pathways")
    print("4. Interaction groups 5. Feedback-loop check 0. Back")

    choice = read_int("Choice: ", 0, 5)

    if choice in (1, 2, 3):
        source = input("Source protein: ").strip().upper()

        if source not in graph:
            print("Error: unknown protein", source)
            return

    if choice == 1:
        print(
            "Recursive DFS:",
            ' -> '.join(dfs_recursive(graph, source))
        )
        print(
            "Iterative DFS:",
            ' -> '.join(dfs_iterative(graph, source))
        )

    elif choice in (2, 3):
        target = input("Target protein: ").strip().upper()

        if target not in graph:
            print("Error: unknown protein", target)
            return

        if choice == 2:
            path = find_pathway(graph, source, target)
            print(
                "Pathway:",
                ' -> '.join(path) if path else "No pathway exists"
            )

        else:
            paths = all_pathways(graph, source, target)

            if not paths:
                print("No pathway exists")

            for number, path in enumerate(paths, 1):
                print(
                    "Pathway",
                    number,
                    ":",
                    ' -> '.join(path)
                )

    elif choice == 4:
        for number, group in enumerate(
            connected_components(graph), 1
        ):
            print(
                "Group",
                number,
                ":",
                ', '.join(group)
            )

    elif choice == 5:
        print("Feedback loop present:", has_cycle(graph))


def run_image_module():
    print("\n--- Module 3: Image-filtering operation cost ---")

    height = read_int("Image height H: ")
    width = read_int("Image width W: ")
    k = read_int("Filter size k (odd): ")

    try:
        out_h, out_w, visits = filter_operation_count(
            height, width, k
        )

        print(
            "Output size:",
            out_h,
            "x",
            out_w,
            "| neighbourhood visits:",
            visits
        )
        print(
            "Formula check:",
            out_h * out_w * k * k
        )

    except ValueError as error:
        print("Error:", error)


def main():
    while True:
        print("\nMAIN MENU")
        print("1. Clinical formula expression tree (BMI, BSA, MAP)")
        print("2. Protein interaction pathways (DFS)")
        print("3. Image-filtering operation analysis")
        print("0. Exit")

        try:
            choice = read_int("Enter choice: ", 0, 3)
        except EOFError:
            break

        if choice == 0:
            print("Exiting program.")
            break

        try:
            [
                run_expression_module,
                run_protein_module,
                run_image_module
            ][choice - 1]()

        except EOFError:
            break


if __name__ == '__main__':
    main()
