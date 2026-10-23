#!/usr/bin/env pymol
'''
PyMOL plugin for US-Flexalign

USAGE:

    usflexalign mobile, target [,args [,exe]]

INSTALLATION

    Install this script as a PyMOL plugin by
    "Plugin" - "Plugin Manager" - "Install New Plugin"

    This plugin depends on the binary executable of US-Flexalign, which must
    be available within a directory specified by PATH, or placed next to this
    plugin file. You can get the PATH value within PyMOL by the following
    command:

    print(os.getenv('PATH'))

    Unlike the "usalign" plugin, this plugin applies a per-hinge-segment
    superposition, so that the flexible alignment result is reproduced in
    PyMOL instead of a single rigid-body transform.
'''
#The host/selection handling in this script follows the usalign plugin by
#Chengxin Zhang (https://github.com/pylelab/USalign), which in turn is partly
#based on the tmalign plugin by Thomas Holder.

from __future__ import print_function

__author__ = 'US-align contributors'
__version__ = '20261023'
__license__ = 'BSD-2-Clause'

from pymol import cmd, CmdException
import subprocess
import tempfile
import shutil
import os
import platform

MATRIX_HEADER = "------ The rotation matrix to rotate Structure_1 to Structure_2 ------"
HINGE_ANNOTATION = "denote different aligned fragment pairs separated by different hinges"

# Same palette order as flexalign_pymol_color_list in
# src/engine/flexalign/flexalign.h, so the plugin and "-pymol" agree on colours.
SEGMENT_COLORS = [
    "blue", "red", "green", "yellow", "magenta", "cyan", "orange", "purple", "brown", "pink",
    "marine", "firebrick", "forest", "gold", "violet", "deepteal", "salmon", "aquamarine", "chocolate", "hotpink",
    "lime", "olive", "slate", "paleblue", "palegreen", "paleyellow", "lightpink", "lightblue", "lightmagenta",
    "dirtyviolet", "dash", "tv_blue", "tv_red", "tv_green",
]
UNALIGNED_COLOR = "grey"
REFERENCE_COLOR = "wheat"


def hinge_char_to_index(hinge_char):
    '''Map a segment marker character to its index, as hinge_char_to_index does
    in src/engine/flexalign/flexalign.h.'''
    if '0' <= hinge_char <= '9':
        return ord(hinge_char) - ord('0')
    if 'a' <= hinge_char <= 'z':
        return ord(hinge_char) - ord('a') + 10
    if 'A' <= hinge_char <= 'Z':
        return ord(hinge_char) - ord('A') + 36
    return 0


def check_executable(exe="USFlexalign"):
    try:
        process = subprocess.Popen('"'+exe+'"', stdout=subprocess.PIPE, shell=True)
        output = process.stdout.read()
        if isinstance(output, bytes):
            output = output.decode('utf-8', 'replace')
        # The US-Flexalign banner spells "Version" with a capital V and the
        # full help spells "version" with a lower-case v, so compare loosely.
        if "version" in output.lower() or "us-flexalign" in output.lower():
            return True
    except OSError:
        print(exe+" not executable")
        return False
    print(exe+" not executable")
    return False


def get_executable_path(exe="USFlexalign"):
    '''Locate the USFlexalign executable.

    Search order: the value of the USFLEXALIGN environment variable, the
    directory holding this plugin, then every directory on PATH.
    '''
    if os.name == "nt" and not exe.lower().endswith(".exe"):
        exe += ".exe"

    candidates = []
    env_exe = os.getenv("USFLEXALIGN", "")
    if env_exe:
        candidates.append(env_exe)
    candidates.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), exe))
    for p in os.getenv("PATH", "").split(os.pathsep):
        if p:
            candidates.append(os.path.join(p, exe))

    for filename in candidates:
        if os.path.isfile(filename) and check_executable(filename):
            return filename

    raise CmdException(
        "Cannot locate '%s'. Place the US-Flexalign executable next to this "
        "plugin, add its directory to PATH, set the USFLEXALIGN environment "
        "variable to its full path, or pass the path as the 4th argument, e.g.\n"
        "    usflexalign mobile, target, , D:\\path\\to\\%s" % (exe, exe))


def _contains_option(args, names):
    tokens = args.split()
    return any(name in tokens for name in names)


def _mirror_enabled(args):
    tokens = args.split()
    for i, token in enumerate(tokens):
        if token == "-mirror" and i + 1 < len(tokens):
            return tokens[i + 1] == "1"
    return False


def _require_loaded(name, role):
    '''Fail early with an actionable message when `name` is not a loaded object.'''
    if not name:
        raise CmdException("usflexalign: missing %s selection" % role)

    count = 0
    try:
        count = cmd.count_atoms("(%s)" % name)
    except Exception:
        count = 0
    if count > 0:
        return

    objects = cmd.get_names("objects")
    loaded = (", ".join(objects) if objects
        else "none - nothing is loaded in this session")

    hint = ""
    if name.lower().endswith((".pdb", ".cif", ".ent", ".pdb1")):
        hint = ("\n'%s' looks like a file path, but this command takes the name of a "
            "structure already loaded in PyMOL. Load it first, e.g.\n"
            "    load %s, my_object" % (name, name))
    raise CmdException(
        "usflexalign: %s selection '%s' contains no atoms.\n"
        "Loaded objects: %s\n"
        "Load the structures first, then refer to them by name:\n"
        "    load <file1.pdb>, mob\n"
        "    load <file2.pdb>, tgt%s" % (role, name, loaded, hint))


def _save_selections(mobile, target, args):
    _require_loaded(mobile, "mobile")
    _require_loaded(target, "target")

    tmpdir = tempfile.mkdtemp(prefix="usflexalign_")
    mobile_filename = os.path.join(tmpdir, "mobile.pdb")
    target_filename = os.path.join(tmpdir, "target.pdb")
    mobile_ca_sele = '(%s) and (not hetatm) and alt +A' % mobile
    target_ca_sele = '(%s) and (not hetatm) and alt +A' % target
    if not "-atom" in args:
        mobile_ca_sele += " and (name CA or name C3')"
        target_ca_sele += " and (name CA or name C3')"

    cmd.save(mobile_filename, mobile_ca_sele)
    cmd.save(target_filename, target_ca_sele)

    # The alignment reported by US-Flexalign covers one representative atom per
    # residue, in the order the residues appear in the saved file. Build the
    # residue list the same way, so that residue index i in seqxA corresponds to
    # residue_keys[i], regardless of how many atoms per residue were saved.
    atoms = []
    cmd.iterate(mobile_ca_sele, "atoms.append((chain, resi))", space={"atoms": atoms})
    residue_keys = []
    seen = set()
    for key in atoms:
        if key in seen:
            continue
        seen.add(key)
        residue_keys.append(key)
    return tmpdir, mobile_filename, target_filename, mobile_ca_sele, residue_keys


def _invoke(exe, mobile_filename, target_filename, args):
    if _contains_option(args, ["-dir", "-dir1", "-dir2", "-dirpair"]):
        raise CmdException("usflexalign: -dir/-dir1/-dir2 cannot be used from the PyMOL plugin")

    if args == '""':
        args = ''
    if len(args) > 2 and args[0] == '"' and args[-1] == '"':
        args = args[1:-1]

    # -outfmt -1 keeps the alignment block (seqxA/seqM/seqyA) and drops the
    # banner; -m - sends the per-segment matrices to stdout. Both are appended
    # last so that they override any value supplied by the user.
    command = ' '.join(['"'+exe+'"', mobile_filename, target_filename, args, '-outfmt -1', '-m -'])
    print(command)

    try:
        process = subprocess.Popen(command, stdout=subprocess.PIPE, shell=True)
        output = process.stdout.read()
    except OSError:
        raise CmdException('Cannot execute "%s", please provide full path to the '
            'US-Flexalign executable' % command)

    if isinstance(output, bytes):
        output = output.decode('utf-8', 'replace')
    return output


def _parse_output(text):
    lines = text.splitlines()

    # ---- rotation matrix blocks: one block per rigid segment ----
    segments = []
    for i, line in enumerate(lines):
        if not line.startswith(MATRIX_HEADER):
            continue
        # the next line is the column caption ("m t[m] u[m][0] ...")
        rows = []
        for k in range(3):
            fields = lines[i + k + 2].split()
            if len(fields) < 5:
                raise CmdException("US-Flexalign: malformed rotation matrix block")
            rows.append([float(x) for x in fields[1:5]])
        segments.append(rows)

    # ---- alignment block: seqxA / seqM / seqyA ----
    annotation = None
    for i, line in enumerate(lines):
        if HINGE_ANNOTATION in line:
            annotation = i
            break
    if annotation is None or annotation + 3 >= len(lines):
        if text.strip():
            # surface the program's own diagnostics (e.g. a rejected -atom value)
            print(text.strip())
        raise CmdException("US-Flexalign: cannot find the alignment block in the program output")

    seqxA = lines[annotation + 1]
    seqM = lines[annotation + 2]

    res2seg = {}
    unaligned = set()
    res_index = -1
    for column in range(min(len(seqxA), len(seqM))):
        if seqxA[column] in '-*':
            continue
        res_index += 1
        hinge_char = seqM[column]
        if hinge_char == ' ':
            unaligned.add(res_index)
            # match find_nearest_hinge_char: search outward for the closest
            # column that carries a segment marker
            hinge_char = '0'
            for offset in range(1, len(seqM)):
                if column - offset >= 0 and seqM[column - offset] != ' ':
                    hinge_char = seqM[column - offset]
                    break
                if column + offset < len(seqM) and seqM[column + offset] != ' ':
                    hinge_char = seqM[column + offset]
                    break
        res2seg[res_index] = hinge_char_to_index(hinge_char)

    if not segments:
        raise CmdException("US-Flexalign: no rotation matrix found in the program output")
    return segments, res2seg, unaligned


def _apply_segment_transforms(object_name, rep_selection, segments, res2seg,
    residue_keys, mirror=False):
    '''Apply one rotation per rigid segment to every atom of the object.

    residue_keys holds the representative atom of each residue in the order the
    program reported them; the selection used to derive it only pins down which
    atom represents a residue, it must not limit which atoms get moved.
    '''
    def seg_of_residue(index):
        segment = res2seg.get(index, 0)
        return segment if segment < len(segments) else 0

    key2seg = {}
    for index, key in enumerate(residue_keys):
        key2seg.setdefault(key, seg_of_residue(index))

    full_sele = 'byobject (%s)' % object_name
    atoms = []
    cmd.iterate_state(1, full_sele, "atoms.append((chain, resi, x, y, z))",
        space={"atoms": atoms})
    if not atoms:
        raise CmdException("usflexalign: object '%s' is empty" % object_name)

    new_coords = []
    for chain, resi, x, y, z in atoms:
        matrix = segments[key2seg.get((chain, resi), 0)]
        if mirror:
            z = -z
        t = [matrix[k][0] for k in range(3)]
        u = [[matrix[k][1 + j] for j in range(3)] for k in range(3)]
        new_coords.append([t[k] + sum(u[k][j] * (x, y, z)[j] for j in range(3))
            for k in range(3)])

    cmd.load_coords(new_coords, full_sele, 1)
    cmd.rebuild()


def _color_by_segment(object_name, segments, res2seg, residue_keys, unaligned,
    color_target=True):
    '''Colour each rigid segment, mirroring the "-pymol" script conventions.

    Segment h gets SEGMENT_COLORS[h]; unaligned residues of the mobile go grey.
    The target is left untouched unless color_target is set, in which case it
    is painted wheat like the reference structure in the generated .pml.
    '''
    seg_to_residues = {}
    for index, key in enumerate(residue_keys):
        segment = res2seg.get(index, 0)
        if segment >= len(segments):
            segment = 0
        if index in unaligned:
            continue
        seg_to_residues.setdefault(segment, []).append(key)

    cmd.color(UNALIGNED_COLOR, object_name)

    made = []
    for segment in sorted(seg_to_residues):
        keys = seg_to_residues[segment]
        selection = _selection_from_keys(object_name, keys)
        name = "hinge_%d_str1" % segment
        cmd.select(name, selection)
        cmd.color(SEGMENT_COLORS[segment % len(SEGMENT_COLORS)], name)
        made.append(name)
    return made


def _selection_from_keys(object_name, keys):
    '''Build a PyMOL selection string for a list of (chain, resi) keys.'''
    by_chain = {}
    for chain, resi in keys:
        by_chain.setdefault(chain, []).append(resi)

    parts = []
    for chain in sorted(by_chain):
        resis = sorted(set(by_chain[chain]))
        if chain:
            parts.append('(%s and chain "%s" and resi %s)'
                % (object_name, chain, '+'.join(resis)))
        else:
            parts.append('(%s and resi %s)' % (object_name, '+'.join(resis)))
    return ' or '.join(parts)


def usflexalign(mobile, target, args='', exe='', color=1):
    '''
USAGE

    usflexalign mobile, target [, args [, exe [, color ]]]

ARGUMENTS

    mobile, target = string: atom selections

    args = string: Extra arguments such as -hinge and -afp

    exe = string: Path to US-Flexalign executable {default: USFlexalign}

    color = 0/1: colour each rigid segment and the reference structure
    {default: 1}

CITATION

    Y Zhu, S Yan, Y Zhang, C Zhang (2026)
    https://github.com/pylelab/USalign
    '''

    # PyMOL passes None for an omitted or blank argument (e.g. "usflexalign a, b, , exe")
    args = args or ''
    exe = exe or ''

    if len(exe) == 0:
        exe = get_executable_path("USFlexalign")
    exe = os.path.abspath(exe) if os.path.isfile(exe) else exe

    tmpdir, mobile_filename, target_filename, mobile_sele, residue_keys = \
        _save_selections(mobile, target, args)
    try:
        output = _invoke(exe, mobile_filename, target_filename, args)
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)

    segments, res2seg, unaligned = _parse_output(output)
    print("US-Flexalign: %d rigid segment(s) detected" % len(segments))

    mapped_residues = (max(res2seg.keys()) + 1) if res2seg else 0
    if len(residue_keys) != mapped_residues:
        print("WARNING! %d representative atoms saved but %d residues mapped; "
            "the superposition may be inexact" % (len(residue_keys), mapped_residues))

    # Every atom of the object moves with its residue's segment, so side chains
    # follow their backbone instead of being left behind.
    objects = cmd.get_object_list("(%s)" % mobile)
    object_name = objects[0] if objects else mobile

    _apply_segment_transforms(object_name, mobile_sele, segments, res2seg,
        residue_keys, mirror=_mirror_enabled(args))

    if int(color):
        target_objects = cmd.get_object_list("(%s)" % target)
        for name in target_objects:
            cmd.color(REFERENCE_COLOR, name)
        made = _color_by_segment(object_name, segments, res2seg, residue_keys, unaligned)
        print("US-Flexalign: coloured %d segment(s): %s"
            % (len(made), ", ".join(made)))


# pymol commands
cmd.extend('usflexalign', usflexalign)
cmd.extend('USFlexalign', usflexalign)

# autocompletion
cmd.auto_arg[0].update({ 'usflexalign': cmd.auto_arg[0]['align'], })
cmd.auto_arg[1].update({ 'usflexalign': cmd.auto_arg[1]['align'], })
cmd.auto_arg[0].update({ 'USFlexalign': cmd.auto_arg[0]['align'], })
cmd.auto_arg[1].update({ 'USFlexalign': cmd.auto_arg[1]['align'], })
