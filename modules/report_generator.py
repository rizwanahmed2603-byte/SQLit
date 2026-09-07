"""
Report generation module for SeqLit.
Renders standalone styled HTML reports and optional PDF reports.
Includes Base64 SVG/Matplotlib visualization generation for offline inclusion in reports.
"""

import os
import json
import base64
from io import BytesIO
from typing import Dict, Any

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    MATPLOTLIB_AVAILABLE = True
except ImportError:
    MATPLOTLIB_AVAILABLE = False

try:
    from weasyprint import HTML
    WEASYPRINT_AVAILABLE = True
except (ImportError, OSError):
    WEASYPRINT_AVAILABLE = False

def generate_static_charts(analysis_data: Dict[str, Any]) -> Dict[str, str]:
    """
    Generates static PNG charts as base64 strings for report inclusion.
    1. Composition chart
    2. BLAST Identity chart
    3. Literature category breakdown
    """
    charts = {}
    if not MATPLOTLIB_AVAILABLE:
        return charts

    try:
        # Chart 1: Base Composition / Frequencies
        freqs = analysis_data.get("validation", {}).get("stats", {}).get("frequencies", {})
        if freqs:
            plt.figure(figsize=(6, 3))
            plt.bar(list(freqs.keys()), list(freqs.values()), color="#206bc4", edgecolor="#18569c")
            plt.title("Sequence Composition (%)", fontsize=10, fontweight="bold")
            plt.ylabel("%")
            plt.ylim(0, 100)
            plt.tight_layout()
            buf = BytesIO()
            plt.savefig(buf, format="png", dpi=120)
            plt.close()
            charts["composition_chart"] = base64.b64encode(buf.getvalue()).decode("utf-8")

        # Chart 2: Literature categories
        articles = analysis_data.get("literature", [])
        if articles:
            cat_counts = {}
            for a in articles:
                c = a.get("category", "General")
                cat_counts[c] = cat_counts.get(c, 0) + 1

            plt.figure(figsize=(5, 3))
            plt.pie(list(cat_counts.values()), labels=list(cat_counts.keys()), autopct="%1.0f%%", startangle=140)
            plt.title("Literature Research Topics", fontsize=10, fontweight="bold")
            plt.tight_layout()
            buf = BytesIO()
            plt.savefig(buf, format="png", dpi=120)
            plt.close()
            charts["literature_chart"] = base64.b64encode(buf.getvalue()).decode("utf-8")

    except Exception as e:
        plt.close("all")

    return charts

def render_html_report(analysis_data: Dict[str, Any], template_str: str) -> str:
    """
    Renders HTML string for the full report using Jinja2.
    """
    from jinja2 import Template
    charts = generate_static_charts(analysis_data)
    template = Template(template_str)
    return template.render(data=analysis_data, charts=charts)

def export_pdf_report(html_content: str, output_path: str) -> bool:
    """
    Converts HTML report to PDF using WeasyPrint if available.
    """
    if not WEASYPRINT_AVAILABLE:
        return False
    try:
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        HTML(string=html_content).write_pdf(output_path)
        return True
    except Exception:
        return False
