#!/usr/bin/env python3
"""
DynoRDG GUI
===========

A themed Tkinter front end for the `dynordg` CLI
(https://github.com/k-meiklejohn/dynordg). It builds the position/event/
probability CSV and the command-line flags for you, runs dynordg as a
subprocess, and renders the resulting figure **inline** in the app (rather
than relying on dynordg's own matplotlib pop-up window).

Requirements
------------
- Python 3.8+ with Tkinter (ships with most standard installs)
- dynordg installed in the same environment this script is run with:
      pip install dynordg
- Pillow (optional, but recommended for a resizable/zoomable PNG preview):
      pip install pillow

Usage
-----
    python dynordg_gui.py

Or, integrated into dynordg's own CLI:

    from dynordg_gui import launch_gui
    # e.g. inside dynordg/__main__.py, when args.gui is set:
    launch_gui()

Design notes
------------
- Output format (PNG or SVG) is chosen explicitly via radio buttons, and the
  correct extension is applied automatically -- you never have to type
  ".svg" yourself. PNG previews inline; SVG can't be rasterised by Tkinter
  without extra dependencies, so it's opened/saved instead.
- dynordg's `--output` flag saves a file without opening a viewer; with no
  `--output` it opens its own interactive matplotlib window. We always pass
  `--output`, defaulting to a temp file when you don't specify a path, so
  the result always ends up embedded in this app.
"""

import os
import re
import sys
import csv
import shutil
import queue
import platform
import threading
import subprocess
import tempfile
import tkinter as tk
from tkinter import ttk, filedialog, messagebox

try:
    from PIL import Image, ImageTk

    HAVE_PIL = True
except ImportError:
    HAVE_PIL = False

EVENT_TYPES = [
    "initiation",
    "termination",
    "40sretention",
    "ires",
    "frameshifting",
    "loadscanning",
    "alldrop",
]

FRAMESHIFTING_RE = re.compile(r"^frameshifting[+-]\d+$")

# ----------------------------------------------------------------- theme --
BG = "#eef1f7"
PANEL_BG = "#ffffff"
ACCENT = "#3a5fe0"
ACCENT_DARK = "#2c48b3"
TEXT = "#1f2430"
MUTED = "#6b7280"
BORDER = "#d7dce6"
LOG_BG = "#12151c"
LOG_FG = "#d7dde5"
LOG_ACCENT = "#8ecbff"


def validate_event_string(event):
    if event in EVENT_TYPES and event != "frameshifting":
        return True
    if FRAMESHIFTING_RE.match(event):
        return True
    return False


def open_path_externally(path):
    """Try to open `path` with the OS default viewer.

    Returns (success: bool, error_message: str | None) so callers can tell
    the user what actually happened instead of failing silently.
    """
    system = platform.system()
    try:
        if system == "Darwin":
            subprocess.run(["open", path], check=True)
        elif system == "Windows":
            os.startfile(path)  # noqa: S606
        else:
            subprocess.run(["xdg-open", path], check=True)
        return True, None
    except FileNotFoundError as exc:
        return False, f"No system opener command was found ({exc})."
    except subprocess.CalledProcessError as exc:
        return False, f"The system opener reported an error (exit code {exc.returncode})."
    except Exception as exc:
        return False, str(exc)


# ------------------------------------------------------------ EventDialog --
class EventDialog(tk.Toplevel):
    """Modal dialog for adding/editing a single (position, event, probability) row."""

    def __init__(self, parent, initial=None):
        super().__init__(parent)
        self.configure(bg=PANEL_BG)
        self.title("Edit event" if initial else "Add event")
        self.resizable(False, False)
        self.result = None
        self.transient(parent)
        self.grab_set()

        pad = {"padx": 8, "pady": 5}
        body = ttk.Frame(self, padding=14, style="Panel.TFrame")
        body.pack(fill="both", expand=True)

        ttk.Label(body, text="Position (nt):", style="Panel.TLabel").grid(row=0, column=0, sticky="e", **pad)
        self.pos_var = tk.StringVar(value=str(initial["position"]) if initial else "")
        ttk.Entry(body, textvariable=self.pos_var, width=14).grid(row=0, column=1, sticky="w", **pad)

        ttk.Label(body, text="Event type:", style="Panel.TLabel").grid(row=1, column=0, sticky="e", **pad)
        base_event = initial["event"] if initial else EVENT_TYPES[0]
        if initial and FRAMESHIFTING_RE.match(initial["event"]):
            base_event = "frameshifting"
        self.event_var = tk.StringVar(value=base_event)
        self.event_combo = ttk.Combobox(body, textvariable=self.event_var, values=EVENT_TYPES, state="readonly", width=16)
        self.event_combo.grid(row=1, column=1, sticky="w", **pad)
        self.event_combo.bind("<<ComboboxSelected>>", self._on_event_change)

        self.fs_frame = ttk.Frame(body, style="Panel.TFrame")
        self.fs_sign_var = tk.StringVar(value="+")
        self.fs_offset_var = tk.StringVar(value="1")
        ttk.Label(self.fs_frame, text="Shift:", style="Panel.TLabel").pack(side="left")
        ttk.Combobox(self.fs_frame, textvariable=self.fs_sign_var, values=["+", "-"], width=3, state="readonly").pack(side="left", padx=(6, 6))
        ttk.Entry(self.fs_frame, textvariable=self.fs_offset_var, width=8).pack(side="left")
        self.fs_frame.grid(row=2, column=0, columnspan=2, sticky="w", padx=8)

        ttk.Label(body, text="Probability:", style="Panel.TLabel").grid(row=3, column=0, sticky="e", **pad)
        self.prob_var = tk.StringVar(value=str(initial["probability"]) if initial else "1")
        ttk.Entry(body, textvariable=self.prob_var, width=14).grid(row=3, column=1, sticky="w", **pad)
        ttk.Label(
            body,
            text="0 < p \u2264 1  (ires / loadscanning may exceed 1)",
            style="Panel.Muted.TLabel",
        ).grid(row=4, column=0, columnspan=2, sticky="w", padx=8)

        if initial and FRAMESHIFTING_RE.match(initial["event"]):
            m = re.match(r"^frameshifting([+-])(\d+)$", initial["event"])
            self.fs_sign_var.set(m.group(1))
            self.fs_offset_var.set(m.group(2))

        btns = ttk.Frame(body, style="Panel.TFrame")
        btns.grid(row=5, column=0, columnspan=2, pady=(14, 0))
        ttk.Button(btns, text="Cancel", command=self.destroy).pack(side="left", padx=6)
        ttk.Button(btns, text="OK", style="Accent.TButton", command=self._on_ok).pack(side="left", padx=6)

        self._on_event_change()
        self.wait_window(self)

    def _on_event_change(self, *_):
        if self.event_var.get() == "frameshifting":
            self.fs_frame.grid()
        else:
            self.fs_frame.grid_remove()

    def _on_ok(self):
        try:
            position = int(self.pos_var.get())
            if position <= 0:
                raise ValueError
        except ValueError:
            messagebox.showerror("Invalid input", "Position must be a positive integer.")
            return

        try:
            probability = float(self.prob_var.get())
        except ValueError:
            messagebox.showerror("Invalid input", "Probability must be a number.")
            return

        if self.event_var.get() == "frameshifting":
            try:
                offset = int(self.fs_offset_var.get())
                if offset <= 0:
                    raise ValueError
            except ValueError:
                messagebox.showerror("Invalid input", "Frameshifting offset must be a positive integer.")
                return
            event = f"frameshifting{self.fs_sign_var.get()}{offset}"
        else:
            event = self.event_var.get()

        if not validate_event_string(event):
            messagebox.showerror("Invalid input", f"'{event}' is not a recognised event type.")
            return

        self.result = {"position": position, "event": event, "probability": probability}
        self.destroy()


# ---------------------------------------------------------------- main app --
class DynoRDGGui(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("DynoRDG")
        self.geometry("1180x780")
        self.minsize(980, 640)
        self.configure(bg=BG)

        self.events = []
        self._log_queue = queue.Queue()
        self._proc = None
        self._last_output_path = None
        self._temp_output_path = None
        self._preview_img_full = None
        self._preview_photo = None
        self._fit_mode = True
        self._zoom = 1.0

        self._setup_style()
        self._build_menu()
        self._build_widgets()
        self.after(120, self._poll_log_queue)

    # --------------------------------------------------------------- style
    def _setup_style(self):
        self.option_add("*Font", ("Segoe UI", 10))
        style = ttk.Style(self)
        try:
            style.theme_use("clam")
        except tk.TclError:
            pass

        style.configure("TFrame", background=BG)
        style.configure("Panel.TFrame", background=PANEL_BG)
        style.configure("TLabel", background=BG, foreground=TEXT)
        style.configure("Panel.TLabel", background=PANEL_BG, foreground=TEXT)
        style.configure("Muted.TLabel", background=BG, foreground=MUTED, font=("Segoe UI", 9))
        style.configure("Panel.Muted.TLabel", background=PANEL_BG, foreground=MUTED, font=("Segoe UI", 9))
        style.configure("Header.TLabel", background=BG, foreground=TEXT, font=("Segoe UI Semibold", 13))
        style.configure("Sub.TLabel", background=BG, foreground=MUTED, font=("Segoe UI", 9))

        style.configure("TLabelframe", background=BG, bordercolor=BORDER, relief="solid")
        style.configure("TLabelframe.Label", background=BG, foreground=TEXT, font=("Segoe UI Semibold", 10))

        style.configure("TButton", padding=6)
        style.configure("Accent.TButton", background=ACCENT, foreground="white", padding=10, font=("Segoe UI Semibold", 11))
        style.map("Accent.TButton", background=[("active", ACCENT_DARK), ("disabled", "#a9b6e8")])

        style.configure("Treeview", background="white", fieldbackground="white", rowheight=24, bordercolor=BORDER)
        style.configure("Treeview.Heading", background="#f0f2f8", foreground=TEXT, font=("Segoe UI Semibold", 9))
        style.map("Treeview", background=[("selected", ACCENT)], foreground=[("selected", "white")])

        style.configure("TEntry", padding=4)
        style.configure("Toolbar.TFrame", background=PANEL_BG)
        style.configure("RunBar.TFrame", background=PANEL_BG)

    # ---------------------------------------------------------------- menu
    def _build_menu(self):
        menubar = tk.Menu(self)

        file_menu = tk.Menu(menubar, tearoff=0)
        file_menu.add_command(label="Load events CSV", command=self._load_csv)
        file_menu.add_command(label="Save events CSV as", command=self._save_csv)
        file_menu.add_separator()
        file_menu.add_command(label="Exit", command=self.destroy)
        menubar.add_cascade(label="File", menu=file_menu)

        run_menu = tk.Menu(menubar, tearoff=0)
        run_menu.add_command(label="Run dynordg", command=self._run)
        run_menu.add_command(label="Show command", command=self._show_command)
        menubar.add_cascade(label="Run", menu=run_menu)

        help_menu = tk.Menu(menubar, tearoff=0)
        help_menu.add_command(label="About", command=self._show_about)
        menubar.add_cascade(label="Help", menu=help_menu)

        self.config(menu=menubar)

    def _show_about(self):
        messagebox.showinfo(
            "About DynoRDG GUI",
            "DynoRDG GUI\n\nA front end for the dynordg CLI.\n"
            "Builds the events CSV and CLI flags, runs dynordg, and shows\n"
            "the resulting figure inline.",
        )

    # ------------------------------------------------------------- layout
    def _build_widgets(self):
        header = ttk.Frame(self, padding=(16, 14, 16, 6))
        header.pack(fill="x")
        ttk.Label(header, text="DynoRDG", style="Header.TLabel").pack(side="left")
        ttk.Label(header, text="  \u2013  dynamic ribosome decision graphs", style="Sub.TLabel").pack(side="left")

        main_pane = ttk.Panedwindow(self, orient="horizontal")
        main_pane.pack(fill="both", expand=True, padx=16, pady=(0, 8))

        left = ttk.Frame(main_pane)
        main_pane.add(left, weight=2)

        right = ttk.Frame(main_pane)
        main_pane.add(right, weight=3)

        self._build_left_panel(left)
        self._build_right_panel(right)

        status_bar = ttk.Frame(self, padding=(16, 4, 16, 10))
        status_bar.pack(fill="x")
        self.status_var = tk.StringVar(value="Ready.")
        ttk.Label(status_bar, textvariable=self.status_var, style="Sub.TLabel").pack(side="left")
        self.progress = ttk.Progressbar(status_bar, mode="indeterminate", length=160)
        self.progress.pack(side="right")

    # -------------------------------------------------------- left panel
    def _build_left_panel(self, parent):
        # --- Run bar: pinned at the bottom of the left column, always visible
        # regardless of how far the panel above is scrolled. Pack this first
        # (side="bottom") so it never gets pushed off-screen.
        run_bar = ttk.Frame(parent, padding=10, style="RunBar.TFrame")
        run_bar.pack(side="bottom", fill="x")
        self.run_btn = ttk.Button(run_bar, text="\u25b6  Run dynordg", style="Accent.TButton", command=self._run)
        self.run_btn.pack(side="left", fill="x", expand=True)
        ttk.Button(run_bar, text="Show command", command=self._show_command).pack(side="left", padx=(8, 0))

        # --- Scrollable content area (events / fasta / parameters / output)
        scroll_container = ttk.Frame(parent)
        scroll_container.pack(side="top", fill="both", expand=True)

        canvas = tk.Canvas(scroll_container, background=BG, highlightthickness=0)
        vscroll = ttk.Scrollbar(scroll_container, orient="vertical", command=canvas.yview)
        canvas.configure(yscrollcommand=vscroll.set)
        canvas.pack(side="left", fill="both", expand=True)
        vscroll.pack(side="right", fill="y")

        inner = ttk.Frame(canvas, padding=(2, 2, 10, 2))
        window_id = canvas.create_window((0, 0), window=inner, anchor="nw")
        inner.bind("<Configure>", lambda e: canvas.configure(scrollregion=canvas.bbox("all")))
        canvas.bind("<Configure>", lambda e: canvas.itemconfig(window_id, width=e.width))
        self._bind_mousewheel(canvas)

        self._build_left_content(inner)

    @staticmethod
    def _bind_mousewheel(canvas):
        def _on_wheel(event):
            if getattr(event, "num", None) == 5 or getattr(event, "delta", 0) < 0:
                canvas.yview_scroll(1, "units")
            else:
                canvas.yview_scroll(-1, "units")

        def _bind(_e=None):
            canvas.bind_all("<MouseWheel>", _on_wheel)
            canvas.bind_all("<Button-4>", _on_wheel)
            canvas.bind_all("<Button-5>", _on_wheel)

        def _unbind(_e=None):
            canvas.unbind_all("<MouseWheel>")
            canvas.unbind_all("<Button-4>")
            canvas.unbind_all("<Button-5>")

        canvas.bind("<Enter>", _bind)
        canvas.bind("<Leave>", _unbind)

    def _build_left_content(self, parent):
        # ---- Events table ----
        events_frame = ttk.LabelFrame(parent, text="Events  (position, event, probability)", padding=8)
        events_frame.pack(fill="x", pady=(0, 10))

        table_row = ttk.Frame(events_frame)
        table_row.pack(fill="x")

        cols = ("position", "event", "probability")
        self.tree = ttk.Treeview(table_row, columns=cols, show="headings", height=8)
        for c, w in zip(cols, (90, 140, 90)):
            self.tree.heading(c, text=c.capitalize())
            self.tree.column(c, width=w, anchor="center")
        self.tree.pack(side="left", fill="x", expand=True)
        self.tree.bind("<Double-1>", lambda e: self._edit_event())

        tree_scroll = ttk.Scrollbar(table_row, orient="vertical", command=self.tree.yview)
        tree_scroll.pack(side="left", fill="y")
        self.tree.configure(yscrollcommand=tree_scroll.set)

        row_btns = ttk.Frame(events_frame)
        row_btns.pack(fill="x", pady=(8, 0))
        ttk.Button(row_btns, text="+ Add", command=self._add_event).pack(side="left")
        ttk.Button(row_btns, text="Edit", command=self._edit_event).pack(side="left", padx=6)
        ttk.Button(row_btns, text="Remove", command=self._remove_event).pack(side="left")
        ttk.Button(row_btns, text="Load CSV", command=self._load_csv).pack(side="right")
        ttk.Button(row_btns, text="Save CSV", command=self._save_csv).pack(side="right", padx=6)

        # ---- FASTA ----
        fasta_frame = ttk.LabelFrame(parent, text="FASTA (optional)", padding=8)
        fasta_frame.pack(fill="x", pady=(0, 10))
        row1 = ttk.Frame(fasta_frame)
        row1.pack(fill="x")
        self.fasta_path_var = tk.StringVar()
        ttk.Entry(row1, textvariable=self.fasta_path_var).pack(side="left", fill="x", expand=True)
        ttk.Button(row1, text="Browse", command=self._browse_fasta).pack(side="left", padx=(6, 0))
        self.guess_var = tk.BooleanVar(value=False)
        ttk.Checkbutton(fasta_frame, text="Guess start sites from Kozak context (-g)", variable=self.guess_var).pack(anchor="w", pady=(6, 0))
        ttk.Label(
            fasta_frame,
            text="Sets transcript length; with -g, experimentally assigns initiation probabilities.",
            style="Muted.TLabel",
            wraplength=380,
        ).pack(anchor="w", pady=(2, 0))

        # ---- Parameters ----
        params_frame = ttk.LabelFrame(parent, text="Simulation parameters", padding=8)
        params_frame.pack(fill="x", pady=(0, 10))
        grid = ttk.Frame(params_frame)
        grid.pack(fill="x")

        def add_row(row, label, var, tip=None, width=12):
            ttk.Label(grid, text=label).grid(row=row, column=0, sticky="e", padx=(0, 8), pady=4)
            ttk.Entry(grid, textvariable=var, width=width).grid(row=row, column=1, sticky="w", pady=4)
            if tip:
                ttk.Label(grid, text=tip, style="Muted.TLabel", wraplength=200).grid(row=row, column=2, sticky="w", padx=(8, 0))

        self.length_var = tk.StringVar()
        add_row(0, "Length (-l):", self.length_var, "default: max position + 10")

        self.logscale_var = tk.StringVar()
        add_row(1, "Log scale (-L):", self.logscale_var, "blank = true-to-scale")

        self.flux_cutoff_var = tk.StringVar()
        add_row(2, "Flux cutoff (-c):", self.flux_cutoff_var, "Redistribute flux below this threshold")

        self.eff_var = tk.StringVar(value="1")
        add_row(3, "5' cap efficiency (-e):", self.eff_var)

        self.trans_decay_var = tk.StringVar()
        add_row(4, "Translation decay (-t):", self.trans_decay_var, "nt half-life")

        self.scan_decay_var = tk.StringVar()
        add_row(5, "Scanning decay (-s):", self.scan_decay_var, "nt half-life")

        self.tc_assoc_var = tk.StringVar()
        add_row(6, "TC re-association (-a):", self.tc_assoc_var, "nt half-life of Ternary complex reassociation during reiniation")

        self.sf_diss_var = tk.StringVar()
        add_row(7, "SF dissociation (-d):", self.sf_diss_var, "nt half-life of scanning factor dissociation after initiation")

        ttk.Separator(grid, orient="horizontal").grid(row=8, column=0, columnspan=3, sticky="ew", pady=8)

        self.init_limit_var = tk.StringVar()
        add_row(9, "Initiation limit (-i):", self.init_limit_var, "experimental / --fasta -g only")





    def _on_format_change(self):
        if self.output_format_var.get() == "png":
            self.output_hint_var.set(
                "The correct .png extension is added automatically. Leave 'Save to' blank "
                "to just preview a temporary file here in the app."
            )
        else:
            self.output_hint_var.set(
                "The correct .svg extension is added automatically. SVG can't be previewed "
                "inline  after running, use 'Open externally' or 'Save output as' below."
            )

    # ------------------------------------------------------- right panel
    def _build_right_panel(self, parent):
        # `parent` is now a plain frame (not a paned window) -- the preview
        # gets the whole right-hand side, since the log lives in its own
        # popup window instead of a docked pane.
        preview_frame = ttk.LabelFrame(parent, text="Output preview", padding=6)
        preview_frame.pack(fill="both", expand=True)

        toolbar = ttk.Frame(preview_frame, style="Toolbar.TFrame")
        toolbar.pack(fill="x", pady=(0, 6))
        ttk.Button(toolbar, text="View log", command=self._show_log_window).pack(side="left")
        ttk.Separator(toolbar, orient="vertical").pack(side="left", fill="y", padx=8)
        ttk.Button(toolbar, text="Fit", command=self._preview_fit).pack(side="left")
        ttk.Button(toolbar, text="100%", command=self._preview_actual_size).pack(side="left", padx=4)
        ttk.Button(toolbar, text="-", width=3, command=lambda: self._preview_zoom_by(1 / 1.25)).pack(side="left")
        ttk.Button(toolbar, text="+", width=3, command=lambda: self._preview_zoom_by(1.25)).pack(side="left", padx=4)
        ttk.Button(toolbar, text="Save output as", command=self._save_output_as).pack(side="right")

        canvas_frame = ttk.Frame(preview_frame)
        canvas_frame.pack(fill="both", expand=True)
        self.canvas = tk.Canvas(canvas_frame, background="#ffffff", highlightthickness=1, highlightbackground=BORDER)
        hbar = ttk.Scrollbar(canvas_frame, orient="horizontal", command=self.canvas.xview)
        vbar = ttk.Scrollbar(canvas_frame, orient="vertical", command=self.canvas.yview)
        self.canvas.configure(xscrollcommand=hbar.set, yscrollcommand=vbar.set)
        self.canvas.grid(row=0, column=0, sticky="nsew")
        vbar.grid(row=0, column=1, sticky="ns")
        hbar.grid(row=1, column=0, sticky="ew")
        canvas_frame.rowconfigure(0, weight=1)
        canvas_frame.columnconfigure(0, weight=1)
        self.canvas.bind("<Configure>", self._on_canvas_resize)
        self._draw_placeholder("Run dynordg to see the graph here")

        if not HAVE_PIL:
            ttk.Label(
                preview_frame,
                text="Tip: install Pillow (pip install pillow) for zoom / fit-to-window support.",
                style="Muted.TLabel",
            ).pack(anchor="w", pady=(4, 0))


                # ---- Output ----
        out_frame = ttk.LabelFrame(parent, text="Output", padding=8)
        out_frame.pack(fill="x", pady=(0, 10))

        fmt_row = ttk.Frame(out_frame)
        fmt_row.pack(fill="x")
        ttk.Label(fmt_row, text="Format:").pack(side="left")
        self.output_format_var = tk.StringVar(value="png")
        ttk.Radiobutton(
            fmt_row, text="PNG  (preview in app)", value="png",
            variable=self.output_format_var, command=self._on_format_change,
        ).pack(side="left", padx=(10, 4))
        ttk.Radiobutton(
            fmt_row, text="SVG  (vector, for editing)", value="svg",
            variable=self.output_format_var, command=self._on_format_change,
        ).pack(side="left", padx=4)

        path_row = ttk.Frame(out_frame)
        path_row.pack(fill="x", pady=(8, 0))
        ttk.Label(path_row, text="Save to:").pack(side="left")
        self.output_var = tk.StringVar()
        ttk.Entry(path_row, textvariable=self.output_var).pack(side="left", fill="x", expand=True, padx=(6, 6))
        ttk.Button(path_row, text="Browse", command=self._browse_output).pack(side="left")

        self.output_hint_var = tk.StringVar()
        ttk.Label(out_frame, textvariable=self.output_hint_var, style="Muted.TLabel", wraplength=380).pack(anchor="w", pady=(6, 0))
        self._on_format_change()

        self._build_log_window()

    # ------------------------------------------------------------ log window
    def _build_log_window(self):
        """Create the log Toplevel once, hidden. `self.log_text` lives here
        for the lifetime of the app and keeps receiving lines in the
        background even while the window is closed/hidden."""
        self.log_window = tk.Toplevel(self)
        self.log_window.title("DynoRDG  Log")
        self.log_window.geometry("640x360")
        self.log_window.configure(bg=BG)
        # Closing the window just hides it -- the log keeps recording.
        self.log_window.protocol("WM_DELETE_WINDOW", self.log_window.withdraw)

        frame = ttk.Frame(self.log_window, padding=8)
        frame.pack(fill="both", expand=True)
        self.log_text = tk.Text(
            frame, wrap="word", state="disabled",
            bg=LOG_BG, fg=LOG_FG, insertbackground=LOG_FG, relief="flat", padx=8, pady=6,
        )
        self.log_text.tag_configure("cmd", foreground=LOG_ACCENT)
        self.log_text.pack(fill="both", expand=True)
        ttk.Button(frame, text="Clear log", command=self._clear_log).pack(anchor="e", pady=(6, 0))

        self.log_window.withdraw()  # start hidden; "View log" brings it up

    def _show_log_window(self):
        if not hasattr(self, "log_window") or not self.log_window.winfo_exists():
            self._build_log_window()
        self.log_window.deiconify()
        self.log_window.lift()
        self.log_window.focus_force()

    # --------------------------------------------------------- preview --
    def _draw_placeholder(self, text):
        self._preview_img_full = None
        self.canvas.delete("all")
        w = max(self.canvas.winfo_width(), 200)
        h = max(self.canvas.winfo_height(), 200)
        self.canvas.create_text(w // 2, h // 2, text=text, fill=MUTED, font=("Segoe UI", 11), justify="center", width=w - 40)
        self.canvas.configure(scrollregion=(0, 0, w, h))

    def _on_canvas_resize(self, _event):
        if self._preview_img_full is None:
            return  # placeholder text doesn't need to be redrawn on every resize
        elif self._fit_mode:
            self._render_preview()

    def _load_preview(self, path):
        if not HAVE_PIL:
            try:
                photo = tk.PhotoImage(file=path)
            except tk.TclError as exc:
                self._log(f"[preview] could not display image natively: {exc}")
                return
            self.canvas.delete("all")
            self.canvas.create_image(0, 0, anchor="nw", image=photo)
            self._preview_photo = photo  # keep reference
            self.canvas.configure(scrollregion=(0, 0, photo.width(), photo.height()))
            return

        try:
            self._preview_img_full = Image.open(path)
            self._preview_img_full.load()
        except Exception as exc:
            self._log(f"[preview] failed to open image: {exc}")
            return
        self._fit_mode = True
        self._render_preview()

    def _render_preview(self):
        if not HAVE_PIL or self._preview_img_full is None:
            return
        img_w, img_h = self._preview_img_full.size
        if self._fit_mode:
            cw = max(self.canvas.winfo_width(), 50)
            ch = max(self.canvas.winfo_height(), 50)
            scale = min(cw / img_w, ch / img_h)
            scale = min(scale, 4.0)
        else:
            scale = self._zoom
        new_w = max(1, int(img_w * scale))
        new_h = max(1, int(img_h * scale))
        resized = self._preview_img_full.resize((new_w, new_h), Image.LANCZOS)
        self._preview_photo = ImageTk.PhotoImage(resized)
        self.canvas.delete("all")
        self.canvas.create_image(0, 0, anchor="nw", image=self._preview_photo)
        self.canvas.configure(scrollregion=(0, 0, new_w, new_h))
        self._zoom = scale

    def _preview_fit(self):
        if self._preview_img_full is None:
            return
        self._fit_mode = True
        self._render_preview()

    def _preview_actual_size(self):
        if self._preview_img_full is None:
            return
        self._fit_mode = False
        self._zoom = 1.0
        self._render_preview()

    def _preview_zoom_by(self, factor):
        if self._preview_img_full is None:
            return
        self._fit_mode = False
        self._zoom = max(0.1, min(self._zoom * factor, 8.0))
        self._render_preview()


    def _save_output_as(self):
        if not self._last_output_path or not os.path.exists(self._last_output_path):
            messagebox.showinfo("No output yet", "Run dynordg first to produce an output file.")
            return
        ext = os.path.splitext(self._last_output_path)[1] or ".png"
        dest = filedialog.asksaveasfilename(defaultextension=ext, filetypes=[("SVG", "*.svg"), ("PNG", "*.png")])
        if not dest:
            return
        try:
            shutil.copy(self._last_output_path, dest)
            self._log(f"Saved output to {dest}")
        except Exception as exc:
            messagebox.showerror("Could not save file", str(exc))

    # ------------------------------------------------------------- events
    def _add_event(self):
        dlg = EventDialog(self)
        if dlg.result:
            self.events.append(dlg.result)
            self._refresh_tree()

    def _edit_event(self):
        sel = self.tree.selection()
        if not sel:
            return
        idx = self.tree.index(sel[0])
        dlg = EventDialog(self, initial=self.events[idx])
        if dlg.result:
            self.events[idx] = dlg.result
            self._refresh_tree()

    def _remove_event(self):
        sel = self.tree.selection()
        if not sel:
            return
        idx = self.tree.index(sel[0])
        del self.events[idx]
        self._refresh_tree()

    def _refresh_tree(self):
        self.events.sort(key=lambda e: e["position"])
        self.tree.delete(*self.tree.get_children())
        for e in self.events:
            self.tree.insert("", "end", values=(e["position"], e["event"], e["probability"]))

    def _load_csv(self):
        path = filedialog.askopenfilename(filetypes=[("CSV files", "*.csv"), ("All files", "*.*")])
        if not path:
            return
        try:
            loaded = []
            with open(path, newline="") as f:
                reader = csv.DictReader(f)
                for row in reader:
                    loaded.append(
                        {
                            "position": int(row["position"]),
                            "event": row["event"].strip(),
                            "probability": float(row["probability"]),
                        }
                    )
        except Exception as exc:
            messagebox.showerror("Could not load CSV", str(exc))
            return
        self.events = loaded
        self._refresh_tree()
        self._log(f"Loaded {len(loaded)} events from {path}")

    def _save_csv(self):
        if not self.events:
            messagebox.showwarning("No events", "There are no events to save yet.")
            return
        path = filedialog.asksaveasfilename(defaultextension=".csv", filetypes=[("CSV files", "*.csv")])
        if not path:
            return
        self._write_csv(path)
        self._log(f"Saved {len(self.events)} events to {path}")

    def _write_csv(self, path):
        with open(path, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["position", "event", "probability"])
            for e in self.events:
                writer.writerow([e["position"], e["event"], e["probability"]])

    # ------------------------------------------------------------- files
    def _browse_fasta(self):
        path = filedialog.askopenfilename(filetypes=[("FASTA files", "*.fasta *.fa *.fna"), ("All files", "*.*")])
        if path:
            self.fasta_path_var.set(path)

    def _browse_output(self):
        fmt = self.output_format_var.get()
        if fmt == "png":
            ext, ftypes = ".png", [("PNG image", "*.png")]
        else:
            ext, ftypes = ".svg", [("SVG image", "*.svg")]
        path = filedialog.asksaveasfilename(defaultextension=ext, filetypes=ftypes)
        if path:
            self.output_var.set(path)

    # ------------------------------------------------------------ command
    def _build_command(self):
        """Returns (cmd_list, tmp_csv_path_or_None, output_path, temp_output_or_None,
        output_format) or raises ValueError."""
        fasta = self.fasta_path_var.get().strip()
        if not self.events and not fasta:
            raise ValueError("Add at least one event, or provide a FASTA file.")

        cmd = [sys.executable, "-m", "dynordg"]
        tmp_csv_path = None

        if self.events:
            tmp_fd, tmp_csv_path = tempfile.mkstemp(suffix=".csv", prefix="dynordg_")
            os.close(tmp_fd)
            self._write_csv(tmp_csv_path)
            cmd.append(tmp_csv_path)

        if fasta:
            cmd += ["--fasta", fasta]
            if self.guess_var.get():
                cmd.append("-g")

        def opt(flag, var):
            val = var.get().strip()
            if val:
                cmd.extend([flag, val])

        opt("-l", self.length_var)
        opt("-L", self.logscale_var)
        opt("-e", self.eff_var)
        opt("-t", self.trans_decay_var)
        opt("-s", self.scan_decay_var)
        opt("-a", self.tc_assoc_var)
        opt("-d", self.sf_diss_var)
        opt("-i", self.init_limit_var)
        opt("-c", self.flux_cutoff_var)

        output_format = self.output_format_var.get()
        ext = ".png" if output_format == "png" else ".svg"
        raw_output = self.output_var.get().strip()
        temp_output = None
        if raw_output:
            # Always normalise the extension to match the chosen format, so a
            # typo or a leftover ".png" while "SVG" is selected can't produce
            # a mismatched file.
            base, _ = os.path.splitext(raw_output)
            output = base + ext
        else:
            fd, temp_output = tempfile.mkstemp(suffix=ext, prefix="dynordg_preview_")
            os.close(fd)
            output = temp_output

        cmd += ["--output", output]
        return cmd, tmp_csv_path, output, temp_output, output_format

    def _show_command(self):
        try:
            cmd, *_ = self._build_command()
        except ValueError as exc:
            messagebox.showerror("Cannot build command", str(exc))
            return
        messagebox.showinfo("Command", " ".join(cmd))

    # ----------------------------------------------------------------run
    def _run(self):
        try:
            cmd, tmp_csv_path, output_path, temp_output, output_format = self._build_command()
        except ValueError as exc:
            messagebox.showerror("Cannot run", str(exc))
            return

        self._log(" ".join(cmd), tag="cmd")
        self.run_btn.config(state="disabled")
        self.status_var.set("Running dynordg")
        self.progress.start(12)

        if self._temp_output_path and os.path.exists(self._temp_output_path):
            try:
                os.remove(self._temp_output_path)
            except OSError:
                pass
        self._temp_output_path = temp_output

        def worker():
            try:
                proc = subprocess.Popen(
                    cmd,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT,
                    text=True,
                    bufsize=1,
                )
                self._proc = proc
                for line in proc.stdout:
                    self._log_queue.put(("line", line.rstrip()))
                proc.wait()
                if proc.returncode == 0:
                    self._log_queue.put(("line", "Finished successfully."))
                    self._log_queue.put(("image", (output_path, output_format)))
                else:
                    self._log_queue.put(("line", f"[process exited with code {proc.returncode}]"))
            except FileNotFoundError:
                self._log_queue.put(
                    ("line", "[error] Could not find dynordg. Install it with: pip install dynordg")
                )
            except Exception as exc:
                self._log_queue.put(("line", f"[error] {exc}"))
            finally:
                if tmp_csv_path:
                    try:
                        os.remove(tmp_csv_path)
                    except OSError:
                        pass
                self._log_queue.put(("done", None))

        threading.Thread(target=worker, daemon=True).start()

    def _poll_log_queue(self):
        try:
            while True:
                kind, payload = self._log_queue.get_nowait()
                if kind == "line":
                    self._log(payload)
                elif kind == "image":
                    output_path, output_format = payload
                    self._last_output_path = output_path
                    if output_format == "png" and os.path.exists(output_path):
                        self._load_preview(output_path)
                    else:
                        self._draw_placeholder(
                            "SVG saved.\nInline preview isn't supported for SVG \n"
                            "use 'Open externally' or 'Save output as' above."
                        )
                        self._log(f"Saved SVG to {output_path}.")
                elif kind == "done":
                    self.run_btn.config(state="normal")
                    self.progress.stop()
                    self.status_var.set("Ready.")
        except queue.Empty:
            pass
        self.after(120, self._poll_log_queue)

    # ------------------------------------------------------------------log
    def _log(self, msg, tag=None):
        self.log_text.config(state="normal")
        if tag:
            self.log_text.insert("end", "$ " + msg + "\n", tag)
        else:
            self.log_text.insert("end", msg + "\n")
        self.log_text.see("end")
        self.log_text.config(state="disabled")

    def _clear_log(self):
        self.log_text.config(state="normal")
        self.log_text.delete("1.0", "end")
        self.log_text.config(state="disabled")


def launch_gui():
    """Entry point for integrating this GUI into another CLI, e.g.:

        if args.gui:
            from dynordg_gui import launch_gui
            launch_gui()
            return
    """
    app = DynoRDGGui()
    app.mainloop()


if __name__ == "__main__":
    launch_gui()