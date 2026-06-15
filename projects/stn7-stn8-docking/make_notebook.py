import nbformat as nbf

nb = nbf.v4.new_notebook()
nb.metadata = {
    "kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
    "language_info": {"name": "python", "version": "3.10.0"},
    "colab": {"name": "STN7_STN8_System_Analysis.ipynb", "provenance": []}
}

cells = []

def md(src): return nbf.v4.new_markdown_cell(src)
def code(src): return nbf.v4.new_code_cell(src)

# ── Title ──
cells.append(md("""# 🖥️ STN7 & STN8: Ubuntu Sistem Analizi
## *Scientific System Performance & Topology Analysis*

---

**Müəllif:** Suleiman Hajizadeh  
**Tarix:** 2026-06-05  
**Sistem:** Ubuntu 26.04 LTS "Resolute" · Kernel 7.0.0-22-generic  
**Mühit:** VMware Virtual Machine · Intel Xeon E5-2620 v4 · 32 vCPU · 31 GB RAM  

---

## 📋 Xülasə (Abstract)

Bu təhlil iki əsas sistem analiz blokunu əhatə edir:

- **STN7** *(System Topology Node 7)* — Sistem topolojisi, NUMA node arxitekturası, CPU nüvə yükü paylanması və yaddaş resurs analizi
- **STN8** *(System Topology Node 8)* — Disk I/O, şəbəkə statistikaları, proses resurs istehlakı və hər CPU nüvəsinin yük profili

> **Elmi əsas:** NUMA (Non-Uniform Memory Access) arxitekturasında yaddaş əlçatanlığı gecikməsi node daxilində 10 ns, nodlar arası 20 ns-dir. Bu fərq yüksək performanslı hesablama iş yüklərinin optimallaşdırılmasında kritik əhəmiyyət daşıyır.
"""))

# ── Install & Import ──
cells.append(md("## 1. Mühit Hazırlığı / Environment Setup"))
cells.append(code("""# Google Colab üçün asılılıqların quraşdırılması
import subprocess, sys

def install(pkg):
    subprocess.check_call([sys.executable, '-m', 'pip', 'install', '-q', pkg])

for p in ['psutil', 'matplotlib', 'numpy', 'seaborn']:
    try:
        __import__(p)
    except ImportError:
        install(p)
        print(f'✅ {p} quraşdırıldı')

print("✅ Bütün kitabxanalar hazırdır")
"""))

cells.append(code("""import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
import matplotlib.ticker as ticker
import numpy as np
import psutil, platform, os, subprocess, time
from datetime import datetime

plt.rcParams.update({
    'font.family': 'DejaVu Sans',
    'axes.unicode_minus': False,
})

# ── Rəng Paleti (Dark Theme) ──────────────────────────────────
C = {
    'bg':     '#0d1117', 'panel':  '#161b22', 'border': '#30363d',
    'blue':   '#58a6ff', 'green':  '#3fb950', 'orange': '#d29922',
    'red':    '#f85149', 'purple': '#bc8cff', 'teal':   '#39d353',
    'text':   '#c9d1d9', 'sub':    '#8b949e',
    'n0':     '#1f6feb', 'n1':     '#388bfd',
}
plt.style.use('dark_background')
print(f"📅 Analiz Tarixi: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
print(f"🐍 Python: {platform.python_version()}")
print(f"📊 Matplotlib: {matplotlib.__version__}")
"""))

# ── Data Collection ──
cells.append(md("""---
## 2. Sistem Məlumatlarının Toplanması
### *Data Acquisition Layer*

Aşağıdakı blok sistemin real-vaxt göstəricilərini `psutil` kitabxanası vasitəsilə toplayır.
"""))

cells.append(code("""def sh(cmd):
    try: return subprocess.check_output(cmd, shell=True, text=True,
                                        stderr=subprocess.DEVNULL).strip()
    except: return "N/A"

# ── CPU ──────────────────────────────────────────────────────
cpu_n       = psutil.cpu_count(logical=True)
cpu_phys    = psutil.cpu_count(logical=False)
cpu_freq    = psutil.cpu_freq()
cpu_pct_all = psutil.cpu_percent(interval=1, percpu=True)
cpu_avg     = np.mean(cpu_pct_all)
load1, load5, load15 = psutil.getloadavg()

# ── Memory ───────────────────────────────────────────────────
mem  = psutil.virtual_memory()
swap = psutil.swap_memory()

# ── Disk ─────────────────────────────────────────────────────
disk = psutil.disk_usage('/')

# ── Network ──────────────────────────────────────────────────
net  = psutil.net_io_counters()

# ── NUMA (hard-coded from numactl --hardware) ─────────────────
numa = {
    0: {'total_mb': 14984, 'free_mb': 11143, 'cpus': list(range(0,16))},
    1: {'total_mb': 16121, 'free_mb': 12310, 'cpus': list(range(16,32))},
}
for nid in numa:
    numa[nid]['used_mb'] = numa[nid]['total_mb'] - numa[nid]['free_mb']
    numa[nid]['pct']     = round(numa[nid]['used_mb'] / numa[nid]['total_mb'] * 100, 1)

# ── System meta ───────────────────────────────────────────────
hostname  = sh("hostname")
os_name   = sh("lsb_release -d").replace("Description:", "").strip()
kernel    = sh("uname -r")
cpu_model = sh("lscpu | grep 'Model name' | cut -d: -f2").strip()
uptime_s  = sh("uptime -p")

print("=" * 55)
print(f"  HOST     : {hostname}")
print(f"  OS       : {os_name}")
print(f"  KERNEL   : {kernel}")
print(f"  CPU      : {cpu_model}")
print(f"  vCPU     : {cpu_n}  |  Ortalama yük: {cpu_avg:.1f}%")
print(f"  RAM      : {mem.total/1024**3:.1f} GB  |  {mem.percent}% istifadə")
print(f"  DISK     : {disk.total/1024**3:.0f} GB  |  {disk.percent}% istifadə")
print(f"  LOAD     : {load1:.2f} / {load5:.2f} / {load15:.2f}")
print(f"  UPTIME   : {uptime_s}")
print("=" * 55)
"""))

# ── STN7 ──
cells.append(md("""---
## 3. STN7 — Sistem Topoloji və Node Baxışı
### *System Topology Node Analysis*

**Elmi izah:**  
STN7 analizi sistemin üst səviyyəli resurs istehlakını, NUMA (Non-Uniform Memory Access) node arxitekturasını
və 32 CPU nüvəsinin yük paylanmasını vizuallaşdırır.

**Əsas göstəricilər:**
- **CPU Utilization Rate:** Ortalama nüvə yükü (%)
- **Memory Pressure:** RAM-ın nə qədərinin istifadə edildiyini göstərir
- **NUMA Imbalance Index:** Node-lar arası yaddaş bərabərsizliyi
- **System Load Average:** 1, 5, 15 dəqiqəlik yük ortalaması (Unix load metric)
"""))

cells.append(code("""fig1, axes = plt.subplots(2, 2, figsize=(16, 11), facecolor=C['bg'])
fig1.suptitle('STN7 — System Topology & Node Overview',
              fontsize=18, fontweight='bold', color=C['text'], y=0.98)

# helper: draw semicircle gauge
def gauge(ax, pct, label, color):
    theta = np.linspace(0, np.pi, 300)
    ax.plot(np.cos(theta), np.sin(theta), color=C['border'], lw=10, solid_capstyle='round')
    filled = np.linspace(0, np.pi * pct / 100, 300)
    ax.plot(np.cos(filled), np.sin(filled), color=color, lw=10, solid_capstyle='round')
    ax.text(0, -0.2, f'{pct:.1f}%', ha='center', fontsize=22,
            fontweight='bold', color=C['text'])
    ax.text(0, -0.52, label, ha='center', fontsize=10, color=C['sub'])
    ax.set_xlim(-1.4,1.4); ax.set_ylim(-0.8,1.25)
    ax.axis('off'); ax.set_facecolor(C['panel'])

# ── A: CPU Gauge ──────────────────────────────────────────────
ax = axes[0,0]
gauge(ax, cpu_avg,
      f'CPU İstifadəsi\\n{cpu_n} nüvə · {cpu_model[:22]}',
      C['green'] if cpu_avg<50 else C['orange'] if cpu_avg<80 else C['red'])
ax.set_title('CPU Utilization Rate', color=C['sub'], fontsize=10)

# ── B: RAM Gauge ──────────────────────────────────────────────
ax = axes[0,1]
gauge(ax, mem.percent,
      f'RAM İstifadəsi\\n{mem.used/1024**3:.1f} GB / {mem.total/1024**3:.0f} GB',
      C['blue'] if mem.percent<50 else C['orange'] if mem.percent<80 else C['red'])
ax.set_title('Memory Utilization Rate', color=C['sub'], fontsize=10)

# ── C: CPU per-core heatmap ──────────────────────────────────
ax = axes[1,0]
ax.set_facecolor(C['panel'])
data = np.array(cpu_pct_all).reshape(2, 16)
im = ax.imshow(data, aspect='auto', cmap='RdYlGn_r', vmin=0, vmax=100)
plt.colorbar(im, ax=ax, label='CPU %', fraction=0.03, pad=0.02)
ax.set_yticks([0,1])
ax.set_yticklabels(['Node 0 (CPU 0–15)', 'Node 1 (CPU 16–31)'],
                   color=C['text'], fontsize=9)
ax.set_xticks(range(16))
ax.set_xticklabels([str(i) for i in range(16)],
                   color=C['sub'], fontsize=8)
for r in range(2):
    for c in range(16):
        v = data[r,c]
        ax.text(c, r, f'{v:.0f}', ha='center', va='center',
                fontsize=7, fontweight='bold',
                color='white' if v>50 else '#111')
ax.set_title('CPU Nüvə Yük Istilik Xəritəsi (32 vCPU)',
             color=C['text'], fontsize=11)

# ── D: NUMA + Load ──────────────────────────────────────────
ax = axes[1,1]
ax.set_facecolor(C['panel'])
# NUMA bars
nodes = ['Node 0\\n(0–15)', 'Node 1\\n(16–31)']
used  = [numa[0]['used_mb']/1024, numa[1]['used_mb']/1024]
free  = [numa[0]['free_mb']/1024, numa[1]['free_mb']/1024]
x = np.arange(2)
ax.bar(x-0.2, used, 0.35, color=C['n0'], label='İstifadə (GB)', alpha=.9)
ax.bar(x+0.2, free, 0.35, color=C['border'], label='Boş (GB)', alpha=.7)
for i,(u,f) in enumerate(zip(used,free)):
    ax.text(i-0.2, u+0.3, f'{u:.1f}', ha='center', color=C['text'], fontsize=9, fontweight='bold')
    ax.text(i+0.2, f+0.3, f'{f:.1f}', ha='center', color=C['sub'], fontsize=9)
# Load avg line (secondary)
ax2 = ax.twinx()
ax2.plot(['Node 0\\n(0–15)', 'Node 1\\n(16–31)'], [load1, load5],
         color=C['orange'], marker='o', lw=2, ms=8, label='Load (1m/5m)')
ax2.set_ylabel('Yük Ortalaması', color=C['orange'], fontsize=9)
ax2.tick_params(colors=C['orange'])
ax.set_xticks(x); ax.set_xticklabels(nodes, color=C['text'])
ax.set_ylabel('GB', color=C['sub']); ax.legend(loc='upper left',
    facecolor=C['panel'], edgecolor=C['border'], labelcolor=C['text'], fontsize=8)
ax.set_title('NUMA Yaddaş Bölgüsü & Yük Ortalaması', color=C['text'], fontsize=11)
ax.tick_params(colors=C['sub']); ax.spines[:].set_color(C['border'])

plt.tight_layout(rect=[0,0,1,0.96])
plt.savefig('STN7_Colab.png', dpi=150, bbox_inches='tight', facecolor=C['bg'])
plt.show()
print("✅ STN7 şəkli saxlandı: STN7_Colab.png")
"""))

# ── STN7 findings ──
cells.append(md("""### 3.1 STN7 Nəticələri / Findings

| Göstərici | Dəyər | Qiymətləndirmə |
|-----------|-------|---------------|
| CPU Utilization | ~21% | ✅ Əla — resurslar bol |
| RAM Pressure | ~10% | ✅ 27 GB sərbəst yaddaş var |
| NUMA Node 0 RAM | 3.8 GB / 14.9 GB | ✅ Balans yaxşıdır |
| NUMA Node 1 RAM | 3.8 GB / 16.1 GB | ✅ Balans yaxşıdır |
| Load Average (1m) | 6.93 | ⚠️ Başlanğıc pikidə — 32 CPU üçün normal |
| Load Average (15m) | 1.36 | ✅ Sistem sabitləşib |

> **Metodoloji qeyd:** Load average dəyəri Unix-də hazır işi gözləyən proseslərin sayını əks etdirir.
> 32 CPU-lu sistemdə 6.93-lük 1-dəqiqəlik yük (~22%) normal start-up davranışıdır.
"""))

# ── STN8 ──
cells.append(md("""---
## 4. STN8 — Disk, Şəbəkə və Proses Analizi
### *Storage, Network & Process Resource Analysis*

**Elmi izah:**  
STN8 analizi sistemin I/O (Input/Output) resurs istehlakını, şəbəkə trafik metrikalarını
və proses səviyyəsində CPU/RAM paylanmasını ölçür.

**Əsas göstəricilər:**
- **Disk Utilization Ratio:** İstifadə olunan disk sahəsinin faizi
- **Network Throughput:** RX/TX bayt/paket sayı
- **Process CPU Affinity:** Hər prosesin nüvə istehlakı
- **Per-Core Load Profile:** 32 nüvənin fərdi yük profili
"""))

cells.append(code("""fig2, axes2 = plt.subplots(2, 2, figsize=(16, 11), facecolor=C['bg'])
fig2.suptitle('STN8 — Storage, Network & Process Analysis',
              fontsize=18, fontweight='bold', color=C['text'], y=0.98)

# ── A: Disk Donut ────────────────────────────────────────────
ax = axes2[0,0]
ax.set_facecolor(C['panel'])
used_g = disk.used/1024**3; free_g = disk.free/1024**3
wedges,_ = ax.pie([used_g, free_g],
                  colors=[C['purple'], C['border']],
                  startangle=90, counterclock=False,
                  wedgeprops=dict(width=0.42, edgecolor=C['bg'], linewidth=2))
ax.text(0, 0.1, f'{used_g:.0f} GB', ha='center', fontsize=18,
        fontweight='bold', color=C['text'])
ax.text(0, -0.18, 'istifadə', ha='center', fontsize=10, color=C['sub'])
ax.text(0, -0.42, f'Cəmi {disk.total/1024**3:.0f} GB  |  Boş {free_g:.0f} GB',
        ha='center', fontsize=8, color=C['sub'])
ax.legend(['İstifadə (%d%%)' % disk.percent, 'Boş (%d%%)' % (100-disk.percent)],
          loc='lower center', facecolor=C['panel'], edgecolor=C['border'],
          labelcolor=C['text'], fontsize=9)
ax.set_title('Disk Utilization Ratio\n/dev/sda2 (SSD)', color=C['text'], fontsize=11)

# ── B: Memory Breakdown ──────────────────────────────────────
ax = axes2[0,1]
ax.set_facecolor(C['panel'])
labels = ['RAM İstifadə', 'Buff/Cache', 'Əlçatan RAM', 'Swap Boş']
vals   = [mem.used/1024**3,
          (mem.buffers+mem.cached)/1024**3,
          mem.available/1024**3,
          swap.total/1024**3]
cols   = [C['blue'], C['teal'], C['border'], C['orange']]
bars   = ax.barh(labels, vals, color=cols, alpha=0.88, height=0.5)
for b, v in zip(bars, vals):
    ax.text(v+0.15, b.get_y()+b.get_height()/2,
            f'{v:.1f} GB', va='center', color=C['text'], fontsize=9, fontweight='bold')
ax.set_xlabel('GB', color=C['sub'])
ax.set_title(f'Memory Breakdown  (Total: {mem.total/1024**3:.0f} GB)',
             color=C['text'], fontsize=11)
ax.tick_params(colors=C['sub']); ax.spines[:].set_color(C['border'])

# ── C: Network I/O ───────────────────────────────────────────
ax = axes2[1,0]
ax.set_facecolor(C['panel'])
net_labels = ['RX Bytes\n(MB)', 'TX Bytes\n(MB)', 'RX Packets\n(K)', 'TX Packets\n(K)']
net_vals   = [net.bytes_recv/1024**2, net.bytes_sent/1024**2,
              net.packets_recv/1000,  net.packets_sent/1000]
raw_vals   = [net.bytes_recv, net.bytes_sent, net.packets_recv, net.packets_sent]
units      = ['MB','MB','paket','paket']
ncols      = [C['green'], C['blue'], C['teal'], C['purple']]
bars_n     = ax.bar(net_labels, net_vals, color=ncols, alpha=0.85, width=0.5)
for b, v, rv, u in zip(bars_n, net_vals, raw_vals, units):
    ax.text(b.get_x()+b.get_width()/2, b.get_height()+max(net_vals)*0.015,
            f'{rv:,}\\n{u}', ha='center', color=C['text'], fontsize=8, fontweight='bold')
ax.set_title('Network Throughput (ens192 · VMware vNIC)',
             color=C['text'], fontsize=11)
ax.tick_params(colors=C['sub']); ax.spines[:].set_color(C['border'])

# ── D: Per-core bar ──────────────────────────────────────────
ax = axes2[1,1]
ax.set_facecolor(C['panel'])
x   = np.arange(cpu_n)
bcl = [C['n0'] if i<16 else C['n1'] for i in range(cpu_n)]
ax.bar(x, cpu_pct_all, color=bcl, alpha=0.85, width=0.75)
ax.axhline(cpu_avg, color=C['orange'], lw=1.8, ls='--',
           label=f'Ortalama: {cpu_avg:.1f}%')
ax.set_xticks(x[::2])
ax.set_xticklabels([str(i) for i in range(0, cpu_n, 2)],
                   color=C['sub'], fontsize=7)
ax.set_ylim(0, 110)
ax.set_xlabel('CPU Nüvə №', color=C['sub'])
ax.set_ylabel('İstifadə %', color=C['sub'])
p0 = mpatches.Patch(color=C['n0'], label='NUMA Node 0 (CPU 0–15)')
p1 = mpatches.Patch(color=C['n1'], label='NUMA Node 1 (CPU 16–31)')
ax.legend(handles=[p0, p1,
          plt.Line2D([0],[0], color=C['orange'], ls='--',
                     label=f'Ort: {cpu_avg:.1f}%')],
          facecolor=C['panel'], edgecolor=C['border'],
          labelcolor=C['text'], fontsize=8)
ax.set_title('Per-Core Load Profile (32 vCPU)', color=C['text'], fontsize=11)
ax.tick_params(colors=C['sub']); ax.spines[:].set_color(C['border'])

plt.tight_layout(rect=[0,0,1,0.96])
plt.savefig('STN8_Colab.png', dpi=150, bbox_inches='tight', facecolor=C['bg'])
plt.show()
print("✅ STN8 şəkli saxlandı: STN8_Colab.png")
"""))

cells.append(md("""### 4.1 STN8 Nəticələri / Findings

| Göstərici | Dəyər | Qiymətləndirmə |
|-----------|-------|---------------|
| Disk Utilization | 43% (40/98 GB) | ✅ Normal — 54 GB boş |
| RAM Buff/Cache | ~5 GB | ✅ OS cache effektivdir |
| Network RX | ~5.1 MB | ✅ Xəta/drop: 0 |
| Network TX | ~9.2 MB | ✅ Sağlam şəbəkə |
| Swap Usage | 0% | ✅ Mükəmməl — RAM kifayət edir |
| Node0/Node1 CPU δ | ~2–5% | ✅ Balans yaxşıdır |

---

## 5. Ümumi Nəticə / Conclusion

```
┌─────────────────────────────────────────────────────┐
│         SİSTEM SƏHHƏTİ DEĞERLENDİRMƏSİ             │
├──────────────┬─────────────┬────────────────────────┤
│  Resurs      │  Status     │  Tövsiyə               │
├──────────────┼─────────────┼────────────────────────┤
│  CPU         │  ✅ Sağlam  │  Normal iş yükü        │
│  RAM         │  ✅ Bol     │  27 GB əlçatan         │
│  Disk        │  ✅ Normal  │  54 GB boş qalır       │
│  Şəbəkə      │  ✅ Sağlam  │  Xəta yoxdur           │
│  Swap        │  ✅ Mükəmməl│  0% istifadə           │
│  NUMA        │  ✅ Balanced│  Node0-Node1 bərabər   │
│  Temperatur  │  ⚠️ Naməlum │  sensors quraşdırın    │
└──────────────┴─────────────┴────────────────────────┘
```

> **Tövsiyə:** `lm-sensors` paketini quraşdıraraq temperatur monitorinqini aktivləşdirin:
> ```bash
> sudo apt install lm-sensors && sudo sensors-detect
> ```
"""))

nb.cells = cells

out = '/home/suleimanhajizadeh/.gemini/antigravity/scratch/colab_analysis/STN7_STN8_Analysis.ipynb'
with open(out, 'w') as f:
    nbf.write(nb, f)

print(f"✅ Notebook yaradıldı: {out}")
