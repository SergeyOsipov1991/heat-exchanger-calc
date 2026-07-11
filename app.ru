import streamlit as st
import math
import io
import os
import matplotlib.pyplot as plt

# Компоненты ReportLab для генерации PDF в оперативной памяти
from reportlab.lib.pagesizes import letter
from reportlab.lib import colors
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, Image
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont

# --- 1. ВСПОМОГАТЕЛЬНЫЕ ФУНКЦИИ И ИНТЕРПОЛЯЦИЯ ПО ГОСТ ---
def interpolate(x, x1, y1, x2, y2):
    """Линейная интерполяция между двумя точками"""
    if x1 == x2:
        return y1
    return y1 + (y2 - y1) * (x - x1) / (x2 - x1)

def get_allowable_stress(material_data, temp):
    """Вычисляет допускаемое напряжение по ГОСТ 34233.1"""
    if temp <= 20:
        return material_data
    points = sorted([t for t in material_data.keys() if t != 20])
    if temp <= points:
        return material_data[points]
    if temp >= points[-1]:
        return material_data[points[-1]]
    for i in range(len(points) - 1):
        t1, t2 = points[i], points[i+1]
        if t1 <= temp <= t2:
            return interpolate(temp, t1, material_data[t1], t2, material_data[t2])
    return material_data[points[-1]]

def get_standard_thickness(calc_thickness):
    """Подбор ближайшей большей стандартной толщины листа по ГОСТ"""
    standard_sizes = [16, 18, 20, 22, 25, 28, 32, 36, 40, 45, 50, 56, 60]
    for size in standard_sizes:
        if size >= calc_thickness:
            return size
    return math.ceil(calc_thickness)

def get_standard_bolt(calc_diameter):
    """Подбор номинального диаметра шпильки по ГОСТ 24705-2004 и ГОСТ 10605-94"""
    bolt_sizes = {
        16: ("M16", 16.0, 2.0),
        20: ("M20", 20.0, 2.5),
        24: ("M24", 24.0, 3.0),
        27: ("M27", 27.0, 3.0),
        30: ("M30", 30.0, 3.5),
        36: ("M36", 36.0, 4.0),
        42: ("M42", 42.0, 4.5),
        48: ("M48", 48.0, 5.0),
        56: ("M56", 56.0, 5.5)
    }
    sorted_keys = sorted(bolt_sizes.keys())
    for size in sorted_keys:
        if size >= calc_diameter:
            return bolt_sizes[size], True
    custom_name = f"M{math.ceil(calc_diameter)} (Превышает М56!)"
    return (custom_name, float(math.ceil(calc_diameter)), 5.5), False

def calculate_full_exchanger(T_rab, selected_material, P_rab, B, A, B_sh, b_p, d_otv, C, z_sh, n_otv, epsilon_friction, A1, B1):
    """Полный цикл прочностного расчета элементов теплообменника"""
    sigma_1_dop_20 = selected_material[20]
    sigma_1_dop_T = get_allowable_stress(selected_material, T_rab)
    
    # Расчет эффективной ширины уплотнения b_э по ГОСТ 34233.4-2017
    b_e = b_p if b_p <= 15.0 else 3.8 * math.sqrt(b_p)
    
    # 1. Расчет проверочного давления
    P_pr = 1.25 * P_rab * (sigma_1_dop_20 / sigma_1_dop_T)
    
    # 2. Коэффициенты плит
    K_m = 0.5 * (B_sh / (B + b_p))
    Y = 1.41 / math.sqrt(1 + (B / A)**2)
    K_o = math.sqrt((1 - (d_otv / B_sh)**3) / (1 - (d_otv / B_sh)))
    
    sigma_1R_dop = sigma_1_dop_T * 1.0918
    
    # 3. Расчет толщин плит
    S1R = K_o * K_m * Y * B * math.sqrt(P_rab / sigma_1R_dop) + C
    S2R = K_o * K_m * Y * B * math.sqrt(P_pr / sigma_1R_dop) + C
    S_H_calc = max(S1R, S2R)
    
    S11R = K_m * Y * B * math.sqrt(P_rab / sigma_1R_dop) + C
    S22R = K_m * Y * B * math.sqrt(P_pr / sigma_1R_dop) + C
    S_pr_calc = max(S11R, S22R)
    
    # 4. Расчет усилий во фланцевом уплотнении
    epsilon = 0.25
    E_n = 6.5
    S_p = 3.0  
    
    q0 = epsilon * (1 + b_e / (2 * S_p)) * E_n
    q_h = 0.8 * 1.2 * P_pr

    F_p = A1 * B1 * P_rab
    F_h = A1 * B1 * P_pr
    F_obzh = 2 * (A1 + B1) * b_p * q0
    F_d = F_obzh + 0.1 * F_p
    F2h = 2 * (A1 + B1) * b_p * q_h

    # 5. Силовой расчет шпилек (до М56)
    sigma_3_dop_T = interpolate(T_rab, 50, 160.0, 200, 130.0)  # Для стали 40Х
    
    P_sh_approx = 3.5
    F_sh_max_approx = max(F_d, F_h + F_p) / 10 
    d_shr_approx = math.sqrt((1.27 * F_sh_max_approx) / (z_sh * sigma_3_dop_T)) * 1000
    
    bolt_info, _ = get_standard_bolt(d_shr_approx)
    bolt_name, d_sh_actual, P_sh_actual = bolt_info
    
    H_sh = 0.5 * P_sh_actual * math.tan(math.radians(60))
    A_sh = (math.pi / 4) * (d_sh_actual - 2 * (17 / 24) * H_sh)**2
    
    F_H = 0.4 * A_sh * z_sh * sigma_3_dop_T
    M_ksh = 0.3 * (max(F_d, F_H) * (d_sh_actual / 1000.0)) / z_sh
    F_sh_z = (M_ksh * z_sh) / (epsilon_friction * (d_sh_actual / 1000.0))
    
    F_sh_p = max(F_d, F_h + F_p, F_sh_z)
    F_sh_h = max(F_d, F2h + F_h, F_sh_z)
    F_sh_max = max(F_sh_p, F_sh_h) / 10

    d_shr = math.sqrt((1.27 * F_sh_max) / (z_sh * sigma_3_dop_T)) * 1000
    
    return S_H_calc, S_pr_calc, P_pr, d_shr, M_ksh, bolt_name, P_sh_actual, d_sh_actual, b_e

def build_pdf_bytes(results_dict, temp_axis, s_h_axis, s_pr_axis):
    """Генерирует PDF-отчет и возвращает его в виде байтового потока."""
    font_path = "arial.ttf"
    if os.path.exists(font_path):
        pdfmetrics.registerFont(TTFont('Arial', font_path))
        font_name = 'Arial'
        font_bold = 'Arial'
    else:
        font_name = 'Helvetica'
        font_bold = 'Helvetica-Bold'

    pdf_buffer = io.BytesIO()
    doc = SimpleDocTemplate(pdf_buffer, pagesize=letter, rightMargin=40, leftMargin=40, topMargin=40, bottomMargin=40)
    styles = getSampleStyleSheet()
    
    title_style = ParagraphStyle('TitleStyle', fontName=font_bold, fontSize=14, leading=17, alignment=1, spaceAfter=12)
    header_style = ParagraphStyle('HeaderStyle', fontName=font_bold, fontSize=11, leading=13, spaceBefore=8, spaceAfter=6)
    normal_style = ParagraphStyle('NormalStyle', fontName=font_name, fontSize=9, leading=13)
    
    story = []
    story.append(Paragraph("ТЕХНИЧЕСКИЙ ОТЧЕТ ПО РАСЧЕТУ НА ПРОЧНОСТЬ И МАССУ", title_style))
    story.append(Paragraph("Разборный пластинчатый теплообменник (ГОСТ 34233 / 33965)", title_style))
    story.append(Spacer(1, 10))
    
    # Таблица 1
    story.append(Paragraph("1. Исходные данные и геометрия уплотнения", header_style))
    data_input = [
        [Paragraph("<b>Параметр конструкции</b>", normal_style), Paragraph("<b>Значение</b>", normal_style)],
        [Paragraph("Материал торцевых плит", normal_style), Paragraph(results_dict["material_name"], normal_style)],
        [Paragraph("Рабочее давление ($P_{раб}$)", normal_style), Paragraph(f"{results_dict['P_rab']} МПа", normal_style)],
        [Paragraph("Расчетная температура ($T$)", normal_style), Paragraph(f"{results_dict['T_rab']} °C", normal_style)],
        [Paragraph("Номинальная ширина прокладки ($b_п$)", normal_style), Paragraph(f"{results_dict['b_p']} мм", normal_style)],
        [Paragraph("Эффективная ширина прокладки по ГОСТ ($b_э$)", normal_style), Paragraph(f"<b>{results_dict['b_e']:.2f} мм</b>", normal_style)],
        [Paragraph("Количество стяжных шпилек ($z_ш$)", normal_style), Paragraph(f"{results_dict['z_sh']} шт", normal_style)],
    ]
    t1 = Table(data_input, colWidths=[240, 240])
    t1.setStyle(TableStyle([
        ('BACKGROUND', (0,0), (1,0), colors.HexColor("#2C3E50")),
        ('TEXTCOLOR', (0,0), (1,0), colors.whitesmoke),
        ('GRID', (0,0), (-1,-1), 0.5, colors.grey),
        ('PADDING', (0,0), (-1,-1), 4),
    ]))
    story.append(t1)
    
    # Таблица 2
    story.append(Paragraph("2. Результаты подбора элементов по ГОСТ сортаментам", header_style))
    data_output = [
        [Paragraph("<b>Наименование компонента</b>", normal_style), 
         Paragraph("<b>Расчетный минимум</b>", normal_style), 
         Paragraph("<b>Подобрано по ГОСТ</b>", normal_style)],
        [Paragraph("Неподвижная плита (толщина)", normal_style), 
         Paragraph(f"{results_dict['S_H_calc']:.2f} мм", normal_style), 
         Paragraph(f"<b>Лист {results_dict['S_H_std']} мм</b>", normal_style)],
        [Paragraph("Прижимная плита (толщина)", normal_style), 
         Paragraph(f"{results_dict['S_pr_calc']:.2f} мм", normal_style), 
         Paragraph(f"<b>Лист {results_dict['S_pr_std']} мм</b>", normal_style)],
        [Paragraph("Стяжные крепежные шпильки (резьба)", normal_style), 
         Paragraph(f"d_шр = {results_dict['d_shr']:.2f} мм", normal_style), 
         Paragraph(f"<b>Резьба {results_dict['bolt_name']}</b>", normal_style)],
    ]
    t2 = Table(data_output, colWidths=[180, 150, 150])
    t2.setStyle(TableStyle([
        ('BACKGROUND', (0,0), (2,0), colors.HexColor("#16A085")),
        ('TEXTCOLOR', (0,0), (2,0), colors.whitesmoke),
        ('GRID', (0,0), (-1,-1), 0.5, colors.grey),
        ('PADDING', (0,0), (-1,-1), 4),
        ('BACKGROUND', (2,1), (2,3), colors.HexColor("#EAEDED")),
    ]))
    story.append(t2)
    
    # Таблица 3
    story.append(Paragraph("3. Весовые характеристики элементов металлоконструкции", header_style))
    data_mass = [
        [Paragraph("<b>Элемент</b>", normal_style), Paragraph("<b>Вес одного элемента, кг</b>", normal_style), Paragraph("<b>Общий вес, кг</b>", normal_style)],
        [Paragraph("Неподвижная плита (с учетом отверстий)", normal_style), Paragraph(f"{results_dict['mass_plate_H']:.1f}", normal_style), Paragraph(f"{results_dict['mass_plate_H']:.1f}", normal_style)],
        [Paragraph("Прижимная плита (сплошная)", normal_style), Paragraph(f"{results_dict['mass_plate_pr']:.1f}", normal_style), Paragraph(f"{results_dict['mass_plate_pr']:.1f}", normal_style)],
        [Paragraph("Комплект стяжного крепежа (шпильки)", normal_style), Paragraph(f"{results_dict['mass_one_bolt']:.2f}", normal_style), Paragraph(f"{results_dict['mass_all_bolts']:.1f}", normal_style)],
        [Paragraph("<b>ИТОГО МЕТАЛЛОКОНСТРУКЦИЯ (без пластин):</b>", normal_style), Paragraph("", normal_style), Paragraph(f"<b>{results_dict['total_mass']:.1f} кг</b>", normal_style)],
    ]
    t3 = Table(data_mass, colWidths=[200, 140, 140])
    t3.setStyle(TableStyle([
        ('BACKGROUND', (0,0), (2,0), colors.HexColor("#D35400")),
        ('TEXTCOLOR', (0,0), (2,0), colors.whitesmoke),
