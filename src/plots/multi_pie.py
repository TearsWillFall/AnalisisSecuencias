import pygal
from pygal.style import Style


# Creates an entry inside a collection.
def create_entry(entry, link_body, color):
    return {
        "value": entry['count'],
        "label": entry.description,
        "xlink": {"href": f"{link_body}{entry.accession}"},
        "color": color
    }


# Creates a list of entries inside a collection.
def create_list(source, min_support, link_body, color):
    result = []

    for index, entry in source.iterrows():
        if (entry['count'] > min_support):
            result.append(create_entry(entry, link_body, color))

    return result


# Plots the multi pie chart and stats.
def plot(go_list, kegg_list, min_support, min_identity, name):
    custom_style = Style(
        opacity='0.8',
        opacity_hover='0.5',
        title_font_size=36,
        tooltip_font_size=10,
        inner_radius=0.75,
        plot_background="rgba(249, 249, 249, 1)"
    )

    multi_pie = pygal.Pie(height=800, tooltip_border_radius=1, style=custom_style)

    go = create_list(go_list, min_support,
                     "https://www.ebi.ac.uk/QuickGO/term/",
                     "rgba(255, 45, 20, .6)")

    kegg = create_list(kegg_list, min_support,
                       "https://www.genome.jp/dbget-bin/www_bget?",
                       "rgba(68, 108, 179, .6)")

    multi_pie.add('GO', go)
    multi_pie.add('KEGG', kegg)

    plot_file_name = f"{name}-out.svg"
    multi_pie.render_to_file(plot_file_name)

    html_file = open(f"{name}.html", "w")
    html_file.write(f"\
<!doctype html>\
<html>\
    <head>\
        <meta charset=\"utf-8\">\
        <title>{name}</title>\
        <style>\
            body {{background-color: #f9f9f9; font-family: Helvetica, Sans-Serif;}}\
            a {{color: blue; text-decoration: none;}}\
        </style>\
    </head>\
    \
    <body>\
        <h1 style=\"text-align: center;\">Functional annotations of {name}</h1>\
        <div style=\"display: flex;\">\
            <object type=\"image/svg+xml\"data=\"{plot_file_name}\" height=\"800\"></object>\
            <div>\
                <div>\
                    <h4>Minimum identity score: {min_identity}</h4>\
                    <h4>Minimum support score: {min_support}</h4>\
                    <div style=\"display: flex;\">\
                        <h2>GO:</h2>\
                        <ul>")

    for go_item in go:
        html_file.write(f"<li><strong>{go_item['value']}x</strong>\
            <a target=\"_blank\" href=\"{go_item['xlink']['href']}\">{go_item['label']}</a></li>")

    html_file.write(f"\
                        </ul>\
                        <h2>KEGG:</h2>\
                        <ul>")

    for kegg_item in kegg:
        html_file.write(f"<li><strong>{kegg_item['value']}x</strong>\
            <a target=\"_blank\" href=\"{kegg_item['xlink']['href']}\">{kegg_item['label']}</a></li>")

    html_file.write(f"\
                        </ul>\
                    </div>\
                </div>\
            </div>\
        </div>\
    </body>\
</html>")
