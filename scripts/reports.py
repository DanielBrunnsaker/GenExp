#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 28 10:47:45 2025

@author: danbru
"""


from pathlib import Path
import tempfile
import pandas as pd
from pdf2image import convert_from_path
from jinja2 import Environment, FileSystemLoader, select_autoescape
from markupsafe import Markup, escape
import markdown
import json

def create_html_report(output_path, sections, template_path,   
                       figure_max_width = "50%"):
    
    """"
    Generate a HTML report with data relevant for the investigation.

    Parameters
    ----------
    output_path : path/string
        Path to the output folder
    sections : list
        List containing relevant informaton needed for the report.
    template_path : path/string
        path to the basic report template used
    figure_max_width : string
        Scaling parameter for the included figures

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
        
        if isinstance(raw_text, (str, Path)) and Path(raw_text).suffix in {'.txt', '.md'}:
            raw_text = Path(raw_text).read_text(encoding='utf-8')
            from_file = True

        if from_file:
            escaped = escape(raw_text)
            html_text = Markup(f'<div style="white-space: pre-wrap; font-family: monospace;">{escaped}</div>')
        else:
            stripped = raw_text.strip()
            if stripped.startswith('<') and stripped.endswith('>'):
                html_text = Markup(raw_text)
            else:
                html_text = Markup(env.filters['markdown'](raw_text))

        # figures
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

        # tables
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



def run_report_generator(experiment_folder):
    
    EXPERIMENT_DIR = Path(experiment_folder)
    
    # Hypothesis text
    file_path = EXPERIMENT_DIR / "hypothesis/selected_hypothesis/hypothesis.txt"

    with open(EXPERIMENT_DIR / "hypothesis/selected_hypothesis/hypothesis_details.json") as f:
        d = json.load(f)
        
    observable = d['observable'].replace('glutamate','glutamic acid')
    
    # rough translation, in case the used words are different

    # Growth data
    df_growth_stats = pd.read_csv(EXPERIMENT_DIR / 'results/growth/tests/AUC_forest.csv')
    df_growth_stats.columns = ['Variable','Beta','SE','CI(2.5%)', 'CI(97.5%)', 'exp(Beta)','Change(%)', 'Pval(empirical)']
    outliers = pd.read_csv(EXPERIMENT_DIR / 'results/growth/outliers.tsv', sep = '\t')
    
    # Target metabolite data
    df_metabolites = pd.read_csv(EXPERIMENT_DIR / 'results/metabolomics/coefficients/supp_stats.csv', sep = '\t')
    df_target = df_metabolites[df_metabolites['Feature'].str.lower().str.contains(observable)]
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
            "text": "Summary of growth experiment. Each subplot corresponds to a the time-series growth data of one well. The shaded areas correspond to the identified log-phase. Plots without any shading (except blanks) are listed as outliers. Pos refers to the intervention compound, and Neg refers to the negative/specificity control.",
            "figures": [(EXPERIMENT_DIR / 'results/plots/growth_curves.pdf', 'Growth curves')],
            "tables": []
        },
        {
            "title": "Metabolomics Data",
            "text": "Overview of metabolomics statistics. Valueare derived from a Mann-Whitney test used on normalized peak areas. FDR adjustment was done with Benjamini-Hochberg.",
            "figures": [],
            "tables": [(df_metabolites, 'Analysis results')]
        },
        {
            "title": "Outliers, growth data",
            "text": "Overview of reported outliers and relevant statistics for growth data. See the methods-section for a summary of how outlier-detection was performed.",
            "figures": [],
            "tables": [(outliers, 'Outlier report, microplate wells.')]
        },
        {
            "title": "Outliers, metabolomics",
            "text": "Overview of reported outliers and relevant statistics for metabolomics data. See the methods-section for a summary of how outlier-detection was performed.",
            "figures": [],
            "tables": [(outliers_metabolomics, 'Outlier report, injections.')]
        }
    ]
    
    create_html_report(
        output_path=EXPERIMENT_DIR / 'report.html',
        sections=sections,
        template_path=EXPERIMENT_DIR / '../../../data/templates/report_template.html',
        figure_max_width='70%'
    )

import argparse

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='create a html report')
    parser.add_argument('--exp_path', required=True)
    args = parser.parse_args()
   
    run_report_generator(args.exp_path)
    




