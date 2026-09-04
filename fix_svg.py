import sys
BS = chr(92)
DQ = chr(34)

def fix_draw_svg_tree(content):
    # Find the function
    start = content.find('static void draw_svg_tree')
    end = content.find('\nbool save_svg')
    if start < 0 or end < 0:
        print("Can't find draw_svg_tree")
        return content
    
    # Build correct function
    lines = []
    def L(s=""):
        lines.append(s)
    
    L("static void draw_svg_tree(")
    L("    ofstream& ofs,")
    L("    const map<string, vector<string>>& children,")
    L("    const map<string, int>& depth,")
    L("    const map<string, double>& ypos,")
    L("    const string& node)")
    L("{")
    L("    const double left_margin = 50;")
    L("    const double x_step = 80;")
    L("    double nx = left_margin + depth.at(node) * x_step;")
    L("    double ny = ypos.at(node);")
    L("    auto it = children.find(node);")
    L("")
    L("    if (it != children.end() && !it->second.empty()) {")
    L("        double first_y = ypos.at(it->second.front());")
    L("        double last_y  = ypos.at(it->second.back());")
    L("")
    L("        if (first_y != last_y)")
    L("            ofs << " + DQ + "<line x1=" + BS + DQ + DQ + " << nx << " + DQ + BS + DQ + " y1=" + BS + DQ + DQ + " << first_y")
    L("                << " + DQ + BS + DQ + " x2=" + BS + DQ + DQ + " << nx << " + DQ + BS + DQ + " y2=" + BS + DQ + DQ + " << last_y")
    L("                << " + DQ + BS + DQ + " stroke=" + BS + DQ + "black" + BS + DQ + " stroke-width=" + BS + DQ + "1" + BS + DQ + "/>" + BS + "n" + DQ + ";")
    
    # Test what one line looks like
    test_line = L  # placeholder
    
    return content  # placeholder

print("Testing output...")
test = DQ + BS + DQ + DQ + " << nx << " + DQ + BS + DQ + " y1=" + BS + DQ + DQ + " << first_y"
print(repr(test))
