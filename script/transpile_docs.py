#!/usr/bin/env python

import os, re

CURRENT_CHAPTER, CURRENT_SECTION, CHAPTER_INDEX, SECTION_INDEX = "", "", 0, 0

with open("tex/main.tex") as tex:

    md_chapter, md_section = None, None

    for line in tex:

        formatted_line = line

        if line.startswith("\chapter"):

            if md_section: md_section.close(); md_section = None

            CURRENT_CHAPTER, CHAPTER_INDEX, SECTION_INDEX = re.search(r"\{(.*)\}\\", line).group(1), CHAPTER_INDEX + 1, 0

            if md_chapter: md_chapter.close()

            md_chapter = open(f"docs/{CURRENT_CHAPTER.lower().replace(' ', '_').replace('--', '_').replace('-', '_')}.md", "w")

            md_chapter.write(f"---\ntitle: {CURRENT_CHAPTER}\nhas_children: true\nlayout: default\nnav_order: {CHAPTER_INDEX}\n---\n\n"); md_chapter.write("{% include mathjax.html %}\n\n")

        if line.startswith("\section"):

            if md_chapter: md_chapter.close(); md_chapter = None

            CURRENT_SECTION, SECTION_INDEX = re.search(r"\{(.*)\}\\", line).group(1), SECTION_INDEX + 1

            if md_section: md_section.close()

            md_section = open(f"docs/{CURRENT_SECTION.lower().replace(' ', '_').replace('--', '_').replace('-', '_')}.md", "w")

            md_section.write(f"---\ntitle: {CURRENT_SECTION}\nparent: {CURRENT_CHAPTER}\nlayout: default\nnav_order: {SECTION_INDEX}\n---\n\n"); md_section.write("{% include mathjax.html %}\n\n")

        if CURRENT_CHAPTER == "" and CURRENT_SECTION == "": continue

        if line.startswith("\chapter") or line.startswith("\section") or line.startswith("\subsection"):

            name = re.search(r"\{(.*)\}\\", line).group(1);

            formatted_line = f"# {name}\n" if line.startswith("\chapter") else f"## {name}\n" if line.startswith("\section") else f"### {name}\n"

        formatted_line = formatted_line.replace("$", "$$").replace("|", "\\\\|")

        formatted_line = formatted_line.replace("\\begin{equation}", "$$\n\\begin{equation}").replace("\\end{equation}", "\\end{equation}\n$$")
        formatted_line = formatted_line.replace("\\begin{align}",    "$$\n\\begin{align}"   ).replace("\\end{align}",    "\\end{align}\n$$"   )

        formatted_line = formatted_line.replace("~\\ref", " \\ref").replace("~\\eqref", " \\eqref")

        if md_chapter: md_chapter.write(f"{formatted_line}")
        if md_section: md_section.write(f"{formatted_line}")
