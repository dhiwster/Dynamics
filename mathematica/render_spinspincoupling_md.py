from __future__ import annotations

from pathlib import Path
import re


SECTION_TITLES = {
    "Single-site eigenvector rotation for the charge-spin blocks",
    "Full symbolic model before lower-charge projection",
    "Rotating-frame transformation and rotating-wave approximation",
    "Symbolic Schrieffer-Wolff generator for the static RWA Hamiltonian",
    "Second-order empty-cavity dispersive Hamiltonian",
    "First-order transformed cavity annihilation operator",
    "Identical-DQD specialization matching the paper",
}


REPLACEMENTS = [
    ("Hd_identical", r"H_{d,identical}"),
    ("J_identical", r"J_{identical}"),
    ("a_eff,vac,identical", r"a_{\mathrm{eff},\mathrm{vac},identical}"),
    ("h_block,-", r"h_{\mathrm{block},-}"),
    ("h_block,+", r"h_{\mathrm{block},+}"),
    ("H0_full", r"H_{0,\mathrm{full}}"),
    ("HI_full", r"H_{I,\mathrm{full}}"),
    ("H0_projected", r"H_{0,\mathrm{proj}}"),
    ("HI_projected", r"H_{I,\mathrm{proj}}"),
    ("HI_rot(t)", r"H_{I,\mathrm{rot}}(t)"),
    ("HI_RWA(t)", r"H_{I,\mathrm{RWA}}(t)"),
    ("HI_RWA", r"H_{I,\mathrm{RWA}}"),
    ("H2_vac", r"H_{2,\mathrm{vac}}"),
    ("Hd_rotating", r"H_{d,\mathrm{rot}}"),
    ("J12", r"J_{12}"),
    ("a_eff,vac", r"a_{\mathrm{eff},\mathrm{vac}}"),
    ("theta_bar", r"\bar{\theta}"),
    ("theta_-", r"\theta_{-}"),
    ("theta_+", r"\theta_{+}"),
    ("E_-", r"E_{-}"),
    ("E_+", r"E_{+}"),
    ("U_-", r"U_{-}"),
    ("U_+", r"U_{+}"),
    ("h_-", r"h_{-}"),
    ("h_+", r"h_{+}"),
    ("omegaQi", r"\omega_{qi}"),
    ("omegaQ1", r"\omega_{q1}"),
    ("omegaQ2", r"\omega_{q2}"),
    ("omegaC", r"\omega_c"),
    ("delta_i", r"\delta_i"),
    ("delta1", r"\delta_1"),
    ("delta2", r"\delta_2"),
    ("Delta1", r"\Delta_1"),
    ("gSigma1", r"g_{\sigma,1}"),
    ("gSigma2", r"g_{\sigma,2}"),
    ("gTau1", r"g_{\tau,1}"),
    ("gTau2", r"g_{\tau,2}"),
    ("gSigma", r"g_{\sigma}"),
    ("gTau", r"g_{\tau}"),
    ("g1", r"g_1"),
    ("g2", r"g_2"),
    ("sp1", r"\sigma_+^{(1)}"),
    ("sm1", r"\sigma_-^{(1)}"),
    ("sz1", r"\sigma_z^{(1)}"),
    ("sp2", r"\sigma_+^{(2)}"),
    ("sm2", r"\sigma_-^{(2)}"),
    ("sz2", r"\sigma_z^{(2)}"),
    ("tau1x", r"\tau_x^{(1)}"),
    ("tau1z", r"\tau_z^{(1)}"),
    ("tau2x", r"\tau_x^{(2)}"),
    ("tau2z", r"\tau_z^{(2)}"),
    ("sigma1x", r"\sigma_x^{(1)}"),
    ("sigma1z", r"\sigma_z^{(1)}"),
    ("sigma2x", r"\sigma_x^{(2)}"),
    ("sigma2z", r"\sigma_z^{(2)}"),
    ("sigma_x", r"\sigma_x"),
    ("ad", r"a^{\dagger}"),
    ("bx", r"b_x"),
    ("Bz", r"B_z"),
    ("tc", r"t_c"),
    ("ETau1", r"E_{\tau,1}"),
    ("ETau2", r"E_{\tau,2}"),
    ("ESigma1", r"E_{\sigma,1}"),
    ("ESigma2", r"E_{\sigma,2}"),
]


def split_top_level_csv(text: str) -> list[str]:
    parts: list[str] = []
    current: list[str] = []
    paren = 0
    brace = 0
    for ch in text:
        if ch == "(":
            paren += 1
        elif ch == ")":
            paren -= 1
        elif ch == "{":
            brace += 1
        elif ch == "}":
            brace -= 1
        elif ch == "," and paren == 0 and brace == 0:
            parts.append("".join(current).strip())
            current = []
            continue
        current.append(ch)
    tail = "".join(current).strip()
    if tail:
        parts.append(tail)
    return parts


def replace_functions(expr: str) -> str:
    def convert_named(expr_in: str, name: str, latex_name: str) -> str:
        out = []
        i = 0
        needle = f"{name}["
        while True:
            j = expr_in.find(needle, i)
            if j == -1:
                out.append(expr_in[i:])
                break
            out.append(expr_in[i:j])
            k = j + len(needle)
            depth = 1
            start = k
            while k < len(expr_in) and depth:
                if expr_in[k] == "[":
                    depth += 1
                elif expr_in[k] == "]":
                    depth -= 1
                k += 1
            inner = expr_in[start : k - 1]
            inner = replace_functions(inner)
            if latex_name == r"\sqrt":
                out.append(rf"{latex_name}{{{inner}}}")
            else:
                out.append(rf"{latex_name}\left({inner}\right)")
            i = k
        return "".join(out)

    expr = convert_named(expr, "ArcTan", r"\arctan")
    expr = convert_named(expr, "Sqrt", r"\sqrt")
    expr = convert_named(expr, "Cos", r"\cos")
    expr = convert_named(expr, "Sin", r"\sin")
    expr = convert_named(expr, "Exp", r"\exp")
    return expr


def replace_parenthesized_superscripts(expr: str) -> str:
    return re.sub(r"\^\(([^()]+)\)", r"^{\1}", expr)


def replace_exponential_notation(expr: str) -> str:
    return re.sub(r"\bE\^\(([^()]+)\)", r"e^{\1}", expr)


def humanize_scalar(expr: str) -> str:
    out = replace_functions(expr.strip())
    for src, dst in REPLACEMENTS:
        out = out.replace(src, dst)
    out = replace_exponential_notation(out)
    out = replace_parenthesized_superscripts(out)
    out = out.replace("^dagger", r"^{\dagger}")
    out = out.replace("**", r"\,")
    out = out.replace("*", " ")
    return out


def humanize_matrix(expr: str) -> str:
    body = expr.strip()
    body = body.removeprefix("{{").removesuffix("}}")
    rows = body.split("}, {")
    parsed_rows = []
    for row in rows:
        clean = row.replace("{", "").replace("}", "")
        parsed_rows.append(split_top_level_csv(clean))
    row_tex = [" & ".join(humanize_scalar(cell) for cell in row) for row in parsed_rows]
    return r"\begin{pmatrix} " + r" \\ ".join(row_tex) + r" \end{pmatrix}"


def humanize_math(expr: str) -> str:
    expr = expr.strip()
    if expr.startswith("{{") and expr.endswith("}}"):
        return humanize_matrix(expr)
    return humanize_scalar(expr)


def render_markdown(lines: list[str]) -> str:
    out: list[str] = ["# Spin-Spin Coupling Derivation", ""]
    for line in lines:
        stripped = line.strip()
        if stripped in SECTION_TITLES:
            out.append(f"## {stripped}")
            out.append("")
            continue
        if not stripped:
            out.append("")
            continue
        if stripped.startswith("NCAlgebra available:"):
            out.extend(["$$", rf"\text{{{stripped}}}", "$$", ""])
            continue
        if stripped.startswith("With paper convention D[c] rho ="):
            out.extend([
                "$$",
                r"\text{With paper convention } \mathcal{D}[c]\,\rho = 2 c \rho c^{\dagger} - c^{\dagger} c\,\rho - \rho\, c^{\dagger} c",
                "$$",
                "",
            ])
            continue
        if " = " in stripped and not stripped.startswith((
            "Following ",
            "RWA drops ",
            "Paper mapping:",
            "The helper ",
            "At the operator level,",
            "Using the same style ",
            "This is the basis-change step ",
            "Using a -> a + [S,a]",
            "Evaluating 1/2 [S, HI_RWA]",
            "Using the standard dispersive SW choice ",
            "Detuning convention:",
            "SW check",
            "With paper convention ",
            "For the static dispersive SW step ",
        )):
            lhs, rhs = stripped.split(" = ", 1)
            out.extend(["$$", f"{humanize_math(lhs)} = {humanize_math(rhs)}", "$$", ""])
            continue
        out.append(humanize_scalar(stripped))
        out.append("")
    return "\n".join(out)


def main() -> None:
    base = Path(__file__).resolve().parent
    txt_path = base / "Spinspincoupling_ncalgebra_paper_derivation_output.txt"
    md_path = base / "Spinspincoupling_ncalgebra_paper_derivation_output.md"
    lines = txt_path.read_text(encoding="utf-8").splitlines()
    md_path.write_text(render_markdown(lines), encoding="utf-8")


if __name__ == "__main__":
    main()
