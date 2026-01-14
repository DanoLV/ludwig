#!/usr/bin/env python3
"""
Simple script to convert markdown to PDF using available tools
"""

import subprocess
import sys
import os

def try_pandoc():
    """Try using pandoc if available"""
    try:
        result = subprocess.run(
            ['pandoc', 'README_velocidad_fluido.md', '-o', 'velocidad_fluido.pdf',
             '--pdf-engine=pdflatex'],
            capture_output=True,
            text=True
        )
        if result.returncode == 0:
            print("PDF created successfully using pandoc")
            return True
        else:
            print(f"Pandoc failed: {result.stderr}")
            return False
    except FileNotFoundError:
        print("Pandoc not found")
        return False

def try_weasyprint():
    """Try using WeasyPrint"""
    try:
        from weasyprint import HTML, CSS
        with open('README_velocidad_fluido.md', 'r') as f:
            content = f.read()

        # Convert markdown to simple HTML
        html_content = markdown_to_html(content)

        HTML(string=html_content).write_pdf('velocidad_fluido.pdf')
        print("PDF created successfully using WeasyPrint")
        return True
    except ImportError:
        print("WeasyPrint not available")
        return False
    except Exception as e:
        print(f"WeasyPrint failed: {e}")
        return False

def markdown_to_html(md_text):
    """Simple markdown to HTML converter"""
    html = """
    <!DOCTYPE html>
    <html>
    <head>
        <meta charset="utf-8">
        <style>
            body {
                font-family: Arial, sans-serif;
                max-width: 800px;
                margin: 40px auto;
                line-height: 1.6;
                padding: 20px;
            }
            h1 { color: #2c3e50; border-bottom: 3px solid #3498db; padding-bottom: 10px; }
            h2 { color: #34495e; margin-top: 30px; border-bottom: 2px solid #95a5a6; padding-bottom: 5px; }
            h3 { color: #7f8c8d; margin-top: 20px; }
            h4 { color: #95a5a6; margin-top: 15px; }
            code {
                background-color: #f4f4f4;
                padding: 2px 6px;
                border-radius: 3px;
                font-family: 'Courier New', monospace;
                font-size: 0.9em;
            }
            pre {
                background-color: #f4f4f4;
                padding: 15px;
                border-radius: 5px;
                overflow-x: auto;
                border-left: 4px solid #3498db;
            }
            pre code {
                background-color: transparent;
                padding: 0;
            }
            table {
                border-collapse: collapse;
                width: 100%;
                margin: 20px 0;
            }
            th, td {
                border: 1px solid #ddd;
                padding: 12px;
                text-align: left;
            }
            th {
                background-color: #3498db;
                color: white;
            }
            td {
                color: #000000;
            }
            tr:nth-child(even) {
                background-color: #f2f2f2;
            }
            strong { color: #2c3e50; }
            hr { border: 0; border-top: 2px solid #ecf0f1; margin: 30px 0; }
        </style>
    </head>
    <body>
    """

    # Simple replacements
    lines = md_text.split('\n')
    in_code_block = False
    in_table = False
    is_first_row = False

    for line in lines:
        # Code blocks
        if line.strip().startswith('```'):
            if in_code_block:
                html += '</code></pre>\n'
                in_code_block = False
            else:
                html += '<pre><code>\n'
                in_code_block = True
            continue

        if in_code_block:
            html += line.replace('<', '&lt;').replace('>', '&gt;') + '\n'
            continue

        # Headers
        if line.startswith('# '):
            html += f'<h1>{line[2:]}</h1>\n'
        elif line.startswith('## '):
            html += f'<h2>{line[3:]}</h2>\n'
        elif line.startswith('### '):
            html += f'<h3>{line[4:]}</h3>\n'
        elif line.startswith('#### '):
            html += f'<h4>{line[5:]}</h4>\n'
        # Horizontal rule
        elif line.strip() == '---':
            html += '<hr>\n'
        # Table
        elif '|' in line and line.strip().startswith('|'):
            if not in_table:
                html += '<table>\n'
                in_table = True
                is_first_row = True

            # Skip separator line
            if '---' in line:
                continue

            cells = [cell.strip() for cell in line.split('|')[1:-1]]

            # First row after table start is the header
            if is_first_row:
                html += '<thead><tr>\n'
                for cell in cells:
                    html += f'<th>{process_inline(cell)}</th>\n'
                html += '</tr></thead>\n<tbody>\n'
                is_first_row = False
            else:
                html += '<tr>\n'
                for cell in cells:
                    html += f'<td>{process_inline(cell)}</td>\n'
                html += '</tr>\n'
        else:
            if in_table and '|' not in line:
                html += '</tbody></table>\n'
                in_table = False
                is_first_row = False

            # Lists
            if line.strip().startswith('- '):
                html += f'<li>{process_inline(line.strip()[2:])}</li>\n'
            # Numbered lists
            elif len(line) > 2 and line[0].isdigit() and line[1:3] == '. ':
                html += f'<li>{process_inline(line.strip()[3:])}</li>\n'
            # Regular paragraphs
            elif line.strip():
                html += f'<p>{process_inline(line)}</p>\n'
            else:
                html += '\n'

    if in_table:
        html += '</tbody></table>\n'

    html += """
    </body>
    </html>
    """
    return html

def process_inline(text):
    """Process inline markdown (bold, italic, code, links)"""
    import re

    # Bold
    text = re.sub(r'\*\*(.+?)\*\*', r'<strong>\1</strong>', text)

    # Code
    text = re.sub(r'`(.+?)`', r'<code>\1</code>', text)

    # Links (simple version)
    text = re.sub(r'\[(.+?)\]\((.+?)\)', r'<a href="\2">\1</a>', text)

    return text

def create_html_fallback():
    """Create an HTML file as fallback"""
    with open('README_velocidad_fluido.md', 'r') as f:
        content = f.read()

    html_content = markdown_to_html(content)

    with open('velocidad_fluido.html', 'w') as f:
        f.write(html_content)

    print("HTML file created: velocidad_fluido.html")
    print("You can open it in a browser and print to PDF (Ctrl+P)")
    return True

if __name__ == '__main__':
    os.chdir('/home/bater/Sim/ludwig/microgel/local')

    # Try different methods
    success = try_pandoc() or try_weasyprint()

    if not success:
        print("\nNo PDF converter found. Creating HTML file instead...")
        create_html_fallback()
        print("\nTo convert HTML to PDF, you can:")
        print("1. Open velocidad_fluido.html in a browser and print to PDF")
        print("2. Install pandoc: sudo apt-get install pandoc texlive-latex-base")
        print("3. Install weasyprint: pip3 install weasyprint")
