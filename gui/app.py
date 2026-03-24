"""Простой русскоязычный GUI для запуска COMTAILS."""
from __future__ import annotations

import io
import os
import threading
from contextlib import redirect_stdout
import tkinter as tk
from tkinter import ttk, messagebox

from simulation import SimulationController
from utils.io_utils import reset_directory
from utils.version import print_version_info, print_citation_info


class ComtailsGUI:
    """Окно управления запуском моделирования COMTAILS."""

    def __init__(self, root: tk.Tk):
        self.root = root
        self.root.title("COMTAILS — русскоязычный интерфейс")
        self.root.geometry("1080x760")

        self.input_dir = tk.StringVar(value="input")
        self.output_dir = tk.StringVar(value="output")
        self.config_file = tk.StringVar(value="TAIL_INPUTS.dat")
        self.dust_profile = tk.StringVar(value="dmdt_vel_power_rmin_rmax.dat")

        self.run_button = None
        self.log_text = None
        self.result_var = tk.StringVar(value="Итог: запуск не выполнялся")
        self.image_label = None
        self.preview_image = None

        self._build_layout()

    def _build_layout(self):
        top = ttk.Frame(self.root, padding=8)
        top.pack(fill=tk.X)

        ttk.Label(top, text="Каталог входных данных:").grid(row=0, column=0, sticky=tk.W, padx=4, pady=4)
        ttk.Entry(top, textvariable=self.input_dir, width=24).grid(row=0, column=1, sticky=tk.W)

        ttk.Label(top, text="Каталог вывода:").grid(row=0, column=2, sticky=tk.W, padx=4, pady=4)
        ttk.Entry(top, textvariable=self.output_dir, width=24).grid(row=0, column=3, sticky=tk.W)

        ttk.Label(top, text="Конфиг:").grid(row=1, column=0, sticky=tk.W, padx=4, pady=4)
        ttk.Entry(top, textvariable=self.config_file, width=24).grid(row=1, column=1, sticky=tk.W)

        ttk.Label(top, text="Профиль пыли:").grid(row=1, column=2, sticky=tk.W, padx=4, pady=4)
        ttk.Entry(top, textvariable=self.dust_profile, width=24).grid(row=1, column=3, sticky=tk.W)

        self.run_button = ttk.Button(top, text="Запустить расчёт", command=self._run_async)
        self.run_button.grid(row=0, column=4, rowspan=2, padx=8, pady=4)

        mid = ttk.PanedWindow(self.root, orient=tk.HORIZONTAL)
        mid.pack(fill=tk.BOTH, expand=True, padx=8, pady=8)

        log_frame = ttk.LabelFrame(mid, text="Журнал вычислений", padding=8)
        self.log_text = tk.Text(log_frame, wrap=tk.WORD, width=72, height=30)
        scroll = ttk.Scrollbar(log_frame, orient=tk.VERTICAL, command=self.log_text.yview)
        self.log_text.configure(yscrollcommand=scroll.set)
        self.log_text.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        scroll.pack(side=tk.RIGHT, fill=tk.Y)
        mid.add(log_frame, weight=3)

        right = ttk.Frame(mid)
        result_frame = ttk.LabelFrame(right, text="Итоговый блок результата", padding=8)
        ttk.Label(result_frame, textvariable=self.result_var, justify=tk.LEFT).pack(anchor=tk.W)
        result_frame.pack(fill=tk.X, pady=(0, 8))

        image_frame = ttk.LabelFrame(right, text="Итоговое изображение", padding=8)
        self.image_label = ttk.Label(image_frame, text="Изображение отсутствует")
        self.image_label.pack(fill=tk.BOTH, expand=True)
        image_frame.pack(fill=tk.BOTH, expand=True)

        mid.add(right, weight=2)

    def _append_log(self, message: str):
        self.log_text.insert(tk.END, message)
        if not message.endswith("\n"):
            self.log_text.insert(tk.END, "\n")
        self.log_text.see(tk.END)

    def _run_async(self):
        self.run_button.configure(state=tk.DISABLED)
        self.result_var.set("Итог: расчёт выполняется...")
        self.log_text.delete("1.0", tk.END)
        thread = threading.Thread(target=self._run_simulation, daemon=True)
        thread.start()

    def _run_simulation(self):
        input_dir = self.input_dir.get().strip()
        output_dir = self.output_dir.get().strip()

        required_files = [
            os.path.join(input_dir, self.config_file.get().strip()),
            os.path.join(input_dir, self.dust_profile.get().strip()),
        ]

        for path in required_files:
            if not os.path.exists(path):
                self.root.after(0, lambda p=path: self._show_error(f"Не найден входной файл: {p}"))
                return

        capture = io.StringIO()

        try:
            with redirect_stdout(capture):
                print_version_info()
                reset_directory(output_dir)
                simulation = SimulationController()
                results = simulation.run(required_files)
                print_version_info()
                print_citation_info()

            self.root.after(0, lambda: self._show_success(capture.getvalue(), results))
        except Exception as exc:  # noqa: BLE001
            log = capture.getvalue()
            self.root.after(0, lambda e=exc, l=log: self._show_exception(e, l))

    def _show_error(self, message: str):
        self._append_log(message)
        messagebox.showerror("Ошибка", message)
        self.run_button.configure(state=tk.NORMAL)

    def _show_exception(self, exc: Exception, log: str):
        if log:
            self._append_log(log)
        self._append_log(f"Критическая ошибка: {exc}")
        self.result_var.set("Итог: выполнение завершилось с ошибкой")
        self.run_button.configure(state=tk.NORMAL)

    def _show_success(self, log: str, results: dict):
        self._append_log(log)

        result_block = (
            "Итог:\n"
            f"• Afrho: {results.get('afrho', 0.0):.4f} м\n"
            f"• Звёздная величина: {results.get('mag', 0.0):.4f}\n"
            f"• Масса выброшенной пыли: {results.get('total_dust_mass', 0.0):.6e} кг\n"
            f"• Изображение: {results.get('image_path') or 'не найдено'}"
        )
        self.result_var.set(result_block)

        self._show_preview(results.get("image_path"))
        self.run_button.configure(state=tk.NORMAL)

    def _show_preview(self, image_path: str | None):
        if not image_path or not os.path.exists(image_path):
            self.image_label.configure(text="Итоговое изображение не найдено", image="")
            self.preview_image = None
            return

        ext = os.path.splitext(image_path)[1].lower()
        try:
            if ext in {".png", ".ppm", ".pgm", ".gif"}:
                self.preview_image = tk.PhotoImage(file=image_path)
                self.image_label.configure(image=self.preview_image, text="")
            else:
                self.image_label.configure(text=f"Файл сформирован: {image_path}", image="")
                self.preview_image = None
        except Exception as exc:  # noqa: BLE001
            self.image_label.configure(text=f"Не удалось показать изображение: {exc}", image="")
            self.preview_image = None


def run_gui():
    root = tk.Tk()
    style = ttk.Style(root)
    if "clam" in style.theme_names():
        style.theme_use("clam")
    app = ComtailsGUI(root)
    root.mainloop()
    return app
