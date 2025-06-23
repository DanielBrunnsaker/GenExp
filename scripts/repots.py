#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 28 10:47:45 2025

@author: danbru
"""


from pathlib import Path
import pandas as pd
from jinja2 import Environment, FileSystemLoader, select_autoescape
import tempfile
from pdf2image import convert_from_path

def create_html_report(
    output_path: Path | str,
    figures: list[Path | tuple[Path | str, str]],
    tables: list[Path | tuple[Path | str, str]],
    template_path: Path | str,    # Full path to the Jinja2 HTML template
    figure_max_width: str = "50%"  # Max width for figures in HTML
):
    """
    Generate an HTML report using a Jinja2 template specified by its file path,
    converting PDFs to PNGs for embedding, and sizing figures.

    :param output_path: Path to write the generated HTML report
    :param figures: List of image paths or (path, caption) tuples
    :param tables: List of CSV paths or (path, caption) tuples
    :param template_path: Path to the Jinja2 HTML template file
    :param figure_max_width: CSS max-width value for <img> tags (e.g., "50%", "400px")
    """
    out_path = Path(output_path)
    tpl_path = Path(template_path)

    # Load Jinja2 template
    env = Environment(
        loader=FileSystemLoader(str(tpl_path.parent)),
        autoescape=select_autoescape([tpl_path.suffix.lstrip('.')])
    )
    template = env.get_template(tpl_path.name)

    # Create temporary directory for converted images
    temp_dir = Path(tempfile.mkdtemp(prefix="report_images_"))

    # Process figures
    fig_data = []
    for item in figures:
        # Unpack path and caption
        if isinstance(item, (list, tuple)) and len(item) == 2:
            img_path, caption = Path(item[0]), item[1]
        else:
            img_path, caption = Path(item), ""

        # Convert PDF to PNG if needed
        if img_path.suffix.lower() == '.pdf':
            pages = convert_from_path(str(img_path), dpi=200, first_page=1, last_page=1)
            png_path = temp_dir / f"{img_path.stem}.png"
            pages[0].save(png_path, "PNG")
            embed_path = png_path.as_posix()
        else:
            embed_path = img_path.as_posix()

        fig_data.append({
            "path": embed_path,
            "caption": caption,
            "max_width": figure_max_width
        })

    # Process tables
    table_data = []
    for item in tables:
        if isinstance(item, (list, tuple)) and len(item) == 2:
            csv_path, caption = Path(item[0]), item[1]
        else:
            csv_path, caption = Path(item), ""
        df = pd.read_csv(csv_path)
        html_tbl = df.to_html(
            classes="report-table",
            index=False,
            border=0,
            justify="left"
        )
        table_data.append({"html": html_tbl, "caption": caption})

    # Render HTML with figure sizing
    rendered = template.render(
        title="My Report",
        figures=fig_data,
        tables=table_data
    )

    # Write output
    out_path.write_text(rendered, encoding='utf-8')



template_path = Path("/Users/danbru/Library/CloudStorage/OneDrive-Chalmers/Desktop/GenExp/data/templates/report_template.html")
template_dir = Path("/Users/danbru/Library/CloudStorage/OneDrive-Chalmers/Desktop/GenExp/data/templates/report_template.html").parent
template_dir.mkdir(parents=True, exist_ok=True)
content = """<!DOCTYPE html>
<html>
<head>…</head>
<body>
  {% for sec in sections %}
    <h2>{{ sec.title }}</h2>
    <p>{{ sec.text }}</p>

    {% for fig in sec.figures %}
      <figure>
        <img src="{{ fig.path }}" style="max-width: {{ fig.max_width }}; height:auto;">
        {% if fig.caption %}<figcaption>{{ fig.caption }}</figcaption>{% endif %}
      </figure>
    {% endfor %}

    {% for tbl in sec.tables %}
      <h3>{{ tbl.caption }}</h3>
      {{ tbl.html | safe }}
    {% endfor %}
  {% endfor %}
</body>
</html>"""
template_path.write_text(content, encoding="utf-8")


fig1_text = """Figure 1. Overview of growth in all of the wells. The color of the border denotes the experimental group. 
The shaded area marks the part of the growth curve used for growth rate calculation. If there is no shaded area, the curve was deemed as an outlier (for non-media-controls.)"""
fig2_text = """Figure 2. Averaged growth curves. Each curve represents the average of the experimental group. All curves were individually smoothed using LOESS. """
fig3_text = """Figure 3. Calculated growth metrics from all of the included wells (Area under growth curve, growth rate and final OD). 
Values were calculated from the processed and smoothed growth curves. Note that growth rate calculations might be not be representative in cases of very slow growth."""

table1_text = """Table 1. Summary statistics calculated from non-parametric bootstrap (5000 resamples). """

figs = [(EXPERIMENT_DIR / "results/plots/growth_curves.pdf", fig1_text),
        (EXPERIMENT_DIR / "results/plots/averaged_curves.pdf", fig2_text),
        EXPERIMENT_DIR / "results/plots/growth_metrics.pdf", fig3_text]

tbls = [(EXPERIMENT_DIR / "results/growth/tests/AUC_forest.csv", "Summary statistics")]


create_html_report(
        output_path=EXPERIMENT_DIR / 'report.html',
        figures=figs,
        tables=tbls,
        template_path=EXPERIMENT_DIR / '../../../data/templates/report_template.html',
        figure_max_width='10px'
    )








from pathlib import Path
import pandas as pd
from jinja2 import Environment, FileSystemLoader, select_autoescape
import tempfile
from pdf2image import convert_from_path


def create_html_report(
    output_path: Path | str,
    sections: list[dict],
    template_path: Path | str,    # Full path to Jinja2 HTML template
    figure_max_width: str = "50%"  # Default max-width for figures
):
    """
    Generate an HTML report with multiple sections. Each section includes free text,
    tables (from CSV or DataFrame), and figures (PDFs converted to PNGs).

    :param output_path: Path to write the generated HTML report
    :param sections: List of dicts with keys:
                     - title: section title (str)
                     - text: free-text content (str)
                     - figures: list of image path or (path, caption)
                     - tables: list of CSV path, DataFrame, or (path/DataFrame, caption)
    :param template_path: Path to the Jinja2 HTML template file
    :param figure_max_width: CSS max-width value for <img> tags (e.g., "50%", "400px")
    """
    out_path = Path(output_path)
    tpl_path = Path(template_path)

    # Load Jinja2 template
    env = Environment(
        loader=FileSystemLoader(str(tpl_path.parent)),
        autoescape=select_autoescape([tpl_path.suffix.lstrip('.')])
    )
    template = env.get_template(tpl_path.name)

    # Temp dir for converting PDFs
    temp_dir = Path(tempfile.mkdtemp(prefix="report_images_"))

    # Prepare section data
    report_sections = []
    for sec in sections:
        sec_title = sec.get("title", "")
        sec_text = sec.get("text", "")

        # Process figures
        sec_figs = []
        for item in sec.get("figures", []):
            if isinstance(item, (list, tuple)) and len(item) == 2:
                img_path, caption = Path(item[0]), item[1]
            else:
                img_path, caption = Path(item), ""
            if img_path.suffix.lower() == '.pdf':
                pages = convert_from_path(str(img_path), dpi=200, first_page=1, last_page=1)
                png_path = temp_dir / f"{img_path.stem}.png"
                pages[0].save(png_path, "PNG")
                embed = png_path.as_posix()
            else:
                embed = img_path.as_posix()
            sec_figs.append({"path": embed, "caption": caption, "max_width": figure_max_width})

        # Process tables (CSV or DataFrame)
        sec_tbls = []
        for item in sec.get("tables", []):
            # Determine DataFrame and caption
            if isinstance(item, (list, tuple)) and len(item) == 2:
                tbl_src, caption = item
            else:
                tbl_src, caption = item, ""
            # Convert to DataFrame if path
            if isinstance(tbl_src, (str, Path)):
                df = pd.read_csv(tbl_src)
            elif isinstance(tbl_src, pd.DataFrame):
                df = tbl_src
            else:
                raise ValueError(f"Unsupported table source: {tbl_src}")
            html_tbl = df.to_html(classes="report-table", index=False, border=0, justify="left")
            sec_tbls.append({"html": html_tbl, "caption": caption})

        report_sections.append({
            "title": sec_title,
            "text": sec_text,
            "figures": sec_figs,
            "tables": sec_tbls
        })

    # Render HTML
    rendered = template.render(sections=report_sections)

    # Write output
    out_path.write_text(rendered, encoding='utf-8')



from pathlib import Path
import tempfile
import pandas as pd
from pdf2image import convert_from_path
from jinja2 import Environment, FileSystemLoader, select_autoescape
from markupsafe import Markup, escape
import markdown

def create_html_report(
    output_path: Path | str,
    sections: list[dict],
    template_path: Path | str,    # Full path to Jinja2 HTML template
    figure_max_width: str = "50%"  # Default max-width for figures
):
    """
    Generate an HTML report with multiple sections. Each section includes structured text,
    tables (from CSV or DataFrame), and figures (PDFs converted to PNGs).

    Text passed as file paths will be embedded exactly as-is, preserving numbering and formatting.

    :param output_path: Path to write the generated HTML report
    :param sections: List of dicts with keys:
                     - title: section title (str)
                     - text: raw text as Markdown/HTML (str) or Path to .txt/.md file
                     - figures: list of image path or (path, caption)
                     - tables: list of CSV path, DataFrame, or (path/DataFrame, caption)
    :param template_path: Path to the Jinja2 HTML template file
    :param figure_max_width: CSS max-width value for <img> tags
    """
    out_path = Path(output_path)
    tpl_path = Path(template_path)

    env = Environment(
        loader=FileSystemLoader(str(tpl_path.parent)),
        autoescape=select_autoescape([tpl_path.suffix.lstrip('.')])
    )
    env.filters['markdown'] = lambda text: markdown.markdown(text or "", extensions=['extra'])
    template = env.get_template(tpl_path.name)

    temp_dir = Path(tempfile.mkdtemp(prefix="report_images_"))
    report_sections = []

    for sec in sections:
        sec_title = sec.get("title", "")
        raw_text = sec.get("text", "")
        from_file = False
        # If text is a file path, read from it and flag
        if isinstance(raw_text, (str, Path)) and Path(raw_text).suffix in {'.txt', '.md'}:
            raw_text = Path(raw_text).read_text(encoding='utf-8')
            from_file = True

        if from_file:
            # Preserve exact formatting (numbers, line breaks) via pre-wrap
            escaped = escape(raw_text)
            html_text = Markup(f'<div style="white-space: pre-wrap; font-family: monospace;">{escaped}</div>')
        else:
            stripped = raw_text.strip()
            if stripped.startswith('<') and stripped.endswith('>'):
                html_text = Markup(raw_text)
            else:
                html_text = Markup(env.filters['markdown'](raw_text))

        # Process figures
        sec_figs = []
        for item in sec.get("figures", []):
            if isinstance(item, (list, tuple)) and len(item) == 2:
                img_path, caption = Path(item[0]), item[1]
            else:
                img_path, caption = Path(item), ""
            if img_path.suffix.lower() == '.pdf':
                pages = convert_from_path(str(img_path), dpi=200, first_page=1, last_page=1)
                png_path = temp_dir / f"{img_path.stem}.png"
                pages[0].save(png_path, "PNG")
                embed = png_path.as_posix()
            else:
                embed = img_path.as_posix()
            sec_figs.append({"path": embed, "caption": caption, "max_width": figure_max_width})

        # Process tables
        sec_tbls = []
        for item in sec.get("tables", []):
            if isinstance(item, (list, tuple)) and len(item) == 2:
                tbl_src, caption = item
            else:
                tbl_src, caption = item, ""
            if isinstance(tbl_src, (str, Path)):
                df = pd.read_csv(tbl_src)
            elif isinstance(tbl_src, pd.DataFrame):
                df = tbl_src
            else:
                raise ValueError(f"Unsupported table source: {tbl_src}")
            html_tbl = df.to_html(classes="report-table", index=False, border=0, justify="left")
            sec_tbls.append({"html": html_tbl, "caption": caption})

        report_sections.append({
            "title": sec_title,
            "text": html_text,
            "figures": sec_figs,
            "tables": sec_tbls
        })

    rendered = template.render(sections=report_sections)
    out_path.write_text(rendered, encoding='utf-8')





if __name__ == '__main__':
    #EXP = Path('.')
    # Define your sections
    
    EXPERIMENT_DIR = Path('/Users/danbru/Library/CloudStorage/OneDrive-Chalmers/Desktop/GenExp/experiments/completed_experiments/aminoadipate_202504291411')
    
    # Hypothesis text
    file_path = EXPERIMENT_DIR / "hypothesis/selected_hypothesis/hypothesis.txt"

    with open(file_path, 'r') as file:
        lines = file.readlines()
        file_content = ''.join(lines)
    
    print(file_content)
    
    # Growth data
    df_growth_stats = pd.read_csv(EXPERIMENT_DIR / 'results/growth/tests/AUC_forest.csv')
    df_growth_stats.columns = ['Variable','Beta','SE','CI(2.5%)', 'CI(97.5%)', 'exp(Beta)','Change(%)', 'Pval(empirical)']
    outliers = pd.read_csv(EXPERIMENT_DIR / 'results/growth/outliers.tsv', sep = '\t')
    
    # Target metabolite data
    df_metabolites = pd.read_csv(EXPERIMENT_DIR / 'results/metabolomics/coefficients/supp_stats.csv', sep = '\t')
    df_target = df_metabolites[df_metabolites['Feature'].str.contains('Aminoadipate')]
    outliers_metabolomics = pd.read_csv(EXPERIMENT_DIR / 'results/metabolomics/processed/outlier_report.tsv', sep = '\t')

    sections = [
        {
            "title": "Hypothesis Summary",
            "text": file_path,
            "figures": [],
            "tables": [(df_growth_stats, 'Summary statistics (GLM)'),
                       (df_target, 'Hypothesis target (Mann-Whitney U)')]
        },
        {
            "title": "Growth Data",
            "text": "Summary of growth experiment. ",
            "figures": [(EXPERIMENT_DIR / 'results/plots/growth_curves.pdf', 'Growth curves')],
            "tables": []
        },
        {
            "title": "Metabolomics Data",
            "text": "Overview of metabolomics results. Statistics are derived from a Mann-Whitney test",
            "figures": [],
            "tables": [(df_metabolites, 'Univariate analysis results (Mann-Whitney U)')]
        },
        {
            "title": "Outliers, growth data",
            "text": "Overview of reported outliers and relevant statistics.",
            "figures": [],
            "tables": [(outliers, 'Outlier report, microplate wells.')]
        },
        {
            "title": "Outliers, metabolomics",
            "text": "Overview of reported outliers and relevant statistics.",
            "figures": [],
            "tables": [(outliers_metabolomics, 'Outlier report, injections.')]
        }
    ]
    
    create_html_report(
        output_path=EXPERIMENT_DIR / 'report.html',
        sections=sections,
        template_path=EXPERIMENT_DIR / '../../../data/templates/report_template.html',
        figure_max_width='50%'
    )




