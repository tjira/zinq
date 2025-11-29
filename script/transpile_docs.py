#!/usr/bin/env python

import os, re

ACRONYMS, CURRENT_CHAPTER, CURRENT_SECTION, CHAPTER_INDEX, SECTION_INDEX, IN_BLOCK, ENUMERATE_INDEX = {}, "", "", 0, 0, False, 0

with open("tex/main.tex") as tex:

    md_chapter, md_section = None, None

    for line in tex:

        formatted_line = line

        if line.startswith(r"\newacronym"):

            acronym = re.findall(r"\{(.*?)\}", line)

            ACRONYMS[acronym[0]] = (acronym[1], acronym[2])

        if line.startswith("\chapter"):

            if md_section: md_section.close(); md_section = None

            CURRENT_CHAPTER, CHAPTER_INDEX, SECTION_INDEX = re.search(r"\{(.*?)\}", line).group(1), CHAPTER_INDEX + 1, 0

            if md_chapter: md_chapter.close()

            md_chapter = open(f"docs/{CURRENT_CHAPTER.lower().replace(' ', '_').replace('--', '_').replace('-', '_')}.md", "w")

            md_chapter.write(f"---\ntitle: {CURRENT_CHAPTER}\nhas_children: true\nlayout: default\nnav_order: {CHAPTER_INDEX}\n---\n\n"); md_chapter.write("{% include mathjax.html %}\n\n")

        if line.startswith("\section"):

            if md_chapter: md_chapter.close(); md_chapter = None

            CURRENT_SECTION, SECTION_INDEX = re.search(r"\{(.*)\}\\", line).group(1), SECTION_INDEX + 1

            if md_section: md_section.close()

            md_section = open(f"docs/{CURRENT_SECTION.lower().replace(' ', '_').replace('--', '_').replace('-', '_')}.md", "w")

            md_section.write(f"---\ntitle: {CURRENT_SECTION}\nparent: {CURRENT_CHAPTER}\nlayout: default\nnav_order: {SECTION_INDEX}\n---\n\n"); md_section.write("{% include mathjax.html %}\n\n")

        if line.startswith(r"\printglossary"): break

        if CURRENT_CHAPTER == "" and CURRENT_SECTION == "": continue

        IN_BLOCK = True  if line.startswith(r"\begin{definition}") else IN_BLOCK
        IN_BLOCK = False if line.startswith(r"\end{definition}"  ) else IN_BLOCK

        if line.startswith(r"\begin{appendices}"): continue
        if line.startswith(r"\end{appendices}"  ): continue

        formatted_line = formatted_line.replace(r"$", r"$$").replace(r"|", r"\\|")

        formatted_line = formatted_line.replace(r"~\ref", r" \ref").replace(r"~\eqref", r" \eqref")

        formatted_line = formatted_line.replace("\\begin{equation}", "\n$$\n\\begin{equation}").replace("\\end{equation}", "\\end{equation}\n$$\n")
        formatted_line = formatted_line.replace("\\begin{align}",    "\n$$\n\\begin{align}"   ).replace("\\end{align}",    "\\end{align}\n$$\n"   )

        if IN_BLOCK: formatted_line = "> " + formatted_line.strip().replace("\n", "\n> ") + "\n"

        formatted_line = "\n{:.definition}\n" if line.startswith(r"\begin{definition}") else formatted_line

        formatted_line = "\n" if line.startswith(r"\end{definition}") else formatted_line

        if line.startswith(r"\subsection"   ): print(line)
        if line.startswith(r"\chapter"      ): formatted_line = "# "    + re.search(r"\{(.*)\}\\", line).group(1) + "\n"
        if line.startswith(r"\section"      ): formatted_line = "## "   + re.search(r"\{(.*)\}\\", line).group(1) + "\n"
        if line.startswith(r"\subsection"   ): formatted_line = "### "  + re.search(r"\{(.*)\}\\", line).group(1) + "\n"
        if line.startswith(r"\subsubsection"): formatted_line = "#### " + re.search(r"\{(.*)\}\\", line).group(1) + "\n"

        for acronym, (short, long) in ACRONYMS.items():
            formatted_line = formatted_line.replace(r"\acrshort{" + acronym + r"}", short).replace(r"\acrfull{" + acronym + r"}", f"{long} ({short})")

        if line.startswith(r"\begin{enumerate}") or line.startswith(r"\end{enumerate}"):
            ENUMERATE_INDEX = 0; continue

        if line.startswith(r"\item"):
            ENUMERATE_INDEX += 1; formatted_line = formatted_line.replace(r"\item", f"{ENUMERATE_INDEX}.")

        if md_chapter: md_chapter.write(f"{formatted_line}")
        if md_section: md_section.write(f"{formatted_line}")
