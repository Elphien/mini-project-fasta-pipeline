"""
=============================================================================
 MINI PROJECT - PIPELINE ANALISIS SEKUENS DNA SEDERHANA
=============================================================================
 Organisme  : Orchidaceae (anggrek) - berbagai spesies genus Cypripedium,
              Paphiopedilum, Phragmipedium, Mexipedium, dll.
 Gen target : 5.8S rRNA gene, ITS1, dan ITS2 (Internal Transcribed Spacer)
 Sumber data: NCBI GenBank (diunduh ulang melalui mirror resmi Biopython,
              https://raw.githubusercontent.com/biopython/biopython/master/
              Doc/examples/ls_orchid.fasta)
              Setiap sekuens memiliki accession number GenBank asli
              (contoh: Z78533.1, Z78532.1, dst.) yang dapat ditelusuri
              langsung di https://www.ncbi.nlm.nih.gov/nuccore/

 Tahapan pipeline (sesuai ketentuan tugas):
   1. Membaca file FASTA
   2. Menyimpan data dalam List
   3. Menghitung frekuensi nukleotida menggunakan Dictionary
   4. Mengurutkan sekuens berdasarkan GC Content
   5. Menampilkan 3 sekuens terbaik (GC Content tertinggi)
   6. Visualisasi grafik hasil berdasarkan nilai GC
   7. Menuliskan hasil ke file CSV
=============================================================================
"""

import csv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# -----------------------------------------------------------------------
# LANGKAH 1 & 2 : MEMBACA FILE FASTA & MENYIMPAN DATA DALAM LIST
# -----------------------------------------------------------------------
def baca_fasta(filepath):
    """
    Membaca file FASTA dan menyimpan setiap record sebagai dictionary
    di dalam sebuah LIST.

    Setiap elemen list berbentuk:
        {
            "id": <ID singkat sekuens>,
            "deskripsi": <deskripsi lengkap header>,
            "sekuens": <string nukleotida ATGC...>
        }

    Parameters
    ----------
    filepath : str
        Path menuju file FASTA (.fasta / .fa)

    Returns
    -------
    list[dict]
        List yang berisi seluruh record sekuens.
    """
    daftar_sekuens = []          # <- LIST utama penyimpanan data
    header = None
    seq_buffer = []

    with open(filepath, "r") as f:
        for baris in f:
            baris = baris.strip()
            if not baris:
                continue
            if baris.startswith(">"):
                # Simpan record sebelumnya (jika ada) sebelum mulai yang baru
                if header is not None:
                    daftar_sekuens.append({
                        "id": header.split()[0][1:],   # buang '>' di depan
                        "deskripsi": header[1:],
                        "sekuens": "".join(seq_buffer).upper()
                    })
                header = baris
                seq_buffer = []
            else:
                seq_buffer.append(baris)

        # Simpan record terakhir setelah loop selesai
        if header is not None:
            daftar_sekuens.append({
                "id": header.split()[0][1:],
                "deskripsi": header[1:],
                "sekuens": "".join(seq_buffer).upper()
            })

    return daftar_sekuens


# -----------------------------------------------------------------------
# LANGKAH 3 : MENGHITUNG FREKUENSI NUKLEOTIDA MENGGUNAKAN DICTIONARY
# -----------------------------------------------------------------------
def hitung_frekuensi_nukleotida(sekuens):
    """
    Menghitung jumlah kemunculan setiap basa nukleotida (A, T, G, C)
    serta basa ambigu lain pada satu sekuens, menggunakan DICTIONARY.

    Parameters
    ----------
    sekuens : str
        String nukleotida.

    Returns
    -------
    dict
        Dictionary frekuensi tiap basa, contoh: {'A': 120, 'T': 98, ...}
    """
    frekuensi = {}   # <- DICTIONARY utama untuk menyimpan frekuensi
    for basa in sekuens:
        frekuensi[basa] = frekuensi.get(basa, 0) + 1
    return frekuensi


def hitung_gc_content(frekuensi, panjang_total):
    """
    Menghitung %GC Content dari dictionary frekuensi nukleotida.

    GC Content (%) = ((jumlah G + jumlah C) / total basa) * 100
    """
    jumlah_g = frekuensi.get("G", 0)
    jumlah_c = frekuensi.get("C", 0)
    if panjang_total == 0:
        return 0.0
    return ((jumlah_g + jumlah_c) / panjang_total) * 100


# -----------------------------------------------------------------------
# PROSES UTAMA PIPELINE
# -----------------------------------------------------------------------
def main():
    input_file = "orchid.fasta"
    output_csv = "hasil_analisis_gc.csv"
    output_plot = "grafik_gc_content.png"

    print("=" * 70)
    print(" PIPELINE ANALISIS SEKUENS DNA - ORCHIDACEAE (DATA NCBI)")
    print("=" * 70)

    # ---- LANGKAH 1 & 2: baca FASTA -> List ----
    data_sekuens = baca_fasta(input_file)
    print(f"\n[1-2] Berhasil membaca {len(data_sekuens)} sekuens dari '{input_file}'")
    print("      Seluruh data disimpan dalam LIST bernama 'data_sekuens'.")

    # ---- LANGKAH 3: hitung frekuensi nukleotida & GC content ----
    print("\n[3]  Menghitung frekuensi nukleotida (Dictionary) & GC Content ...")
    for record in data_sekuens:
        seq = record["sekuens"]
        freq = hitung_frekuensi_nukleotida(seq)         # Dictionary
        gc = hitung_gc_content(freq, len(seq))

        record["panjang"] = len(seq)
        record["frekuensi"] = freq
        record["jumlah_A"] = freq.get("A", 0)
        record["jumlah_T"] = freq.get("T", 0)
        record["jumlah_G"] = freq.get("G", 0)
        record["jumlah_C"] = freq.get("C", 0)
        record["gc_content"] = round(gc, 2)

    # Tampilkan contoh dictionary frekuensi salah satu sekuens
    contoh = data_sekuens[0]
    print(f"      Contoh dictionary frekuensi untuk '{contoh['id']}':")
    print(f"      {contoh['frekuensi']}")

    # ---- LANGKAH 4: urutkan berdasarkan GC Content (descending) ----
    data_terurut = sorted(data_sekuens, key=lambda x: x["gc_content"], reverse=True)
    print("\n[4]  Sekuens telah diurutkan berdasarkan GC Content (tertinggi -> terendah).")

    # ---- LANGKAH 5: tampilkan 3 sekuens terbaik ----
    print("\n[5]  TOP 3 SEKUENS DENGAN GC CONTENT TERTINGGI:")
    print("-" * 70)
    top3 = data_terurut[:3]
    for i, rec in enumerate(top3, start=1):
        print(f"  Peringkat {i}")
        print(f"    ID            : {rec['id']}")
        print(f"    Deskripsi     : {rec['deskripsi']}")
        print(f"    Panjang       : {rec['panjang']} bp")
        print(f"    GC Content    : {rec['gc_content']}%")
        print(f"    Komposisi     : A={rec['jumlah_A']}  T={rec['jumlah_T']}  "
              f"G={rec['jumlah_G']}  C={rec['jumlah_C']}")
        print("-" * 70)

    # ---- LANGKAH 6: visualisasi grafik GC Content ----
    print("\n[6]  Membuat visualisasi grafik GC Content ...")
    buat_grafik(data_terurut, top3, output_plot)
    print(f"      Grafik tersimpan sebagai '{output_plot}'")

    # ---- LANGKAH 7: tulis hasil ke file CSV ----
    print("\n[7]  Menuliskan seluruh hasil analisis ke file CSV ...")
    tulis_csv(data_terurut, output_csv)
    print(f"      Data tersimpan sebagai '{output_csv}'")

    print("\n" + "=" * 70)
    print(" PIPELINE SELESAI DIJALANKAN")
    print("=" * 70)

    return data_terurut, top3


def buat_grafik(data_terurut, top3, output_path):
    """
    Membuat figure dengan 2 panel:
      Panel kiri  : bar chart seluruh sekuens (diurutkan), top-3 disorot
      Panel kanan : bar chart horizontal khusus 3 sekuens terbaik (rinci)
    """
    ids = [rec["id"].split("|")[-1] for rec in data_terurut]   # pakai kode pendek, misal Z78533.1|CIZ78533
    gc_values = [rec["gc_content"] for rec in data_terurut]
    top3_ids = {rec["id"] for rec in top3}
    colors = ["#2E8B57" if rec["id"] in top3_ids else "#A9D6B5" for rec in data_terurut]
    rata2 = sum(gc_values) / len(gc_values)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(17, 7), gridspec_kw={"width_ratios": [2.3, 1]})

    # ---- Panel 1: seluruh sekuens ----
    ax1.bar(range(len(ids)), gc_values, color=colors, edgecolor="white", linewidth=0.3)
    ax1.set_title("GC Content Seluruh Sekuens Orchidaceae (n=94)\nData: NCBI GenBank",
                   fontsize=12, fontweight="bold")
    ax1.set_xlabel("Sekuens, diurutkan menurun berdasarkan GC Content", fontsize=10)
    ax1.set_ylabel("GC Content (%)", fontsize=10)
    ax1.set_xticks([])
    ax1.axhline(y=rata2, color="gray", linestyle="--", linewidth=1,
                label=f"Rata-rata = {rata2:.2f}%")
    ax1.set_ylim(0, max(gc_values) + 5)
    ax1.legend(loc="upper right", fontsize=9)

    # Tandai 3 batang teratas dengan garis penunjuk kecil (tanpa numpuk teks)
    for idx, rec in enumerate(top3):
        pos = [i for i, r in enumerate(data_terurut) if r["id"] == rec["id"]][0]
        ax1.annotate(f"#{idx+1}", xy=(pos, rec["gc_content"]),
                     xytext=(pos, rec["gc_content"] + 2.2 + idx * 1.4),
                     ha="center", fontsize=9, fontweight="bold", color="#1B5E3A",
                     arrowprops=dict(arrowstyle="-", color="#1B5E3A", lw=0.8))

    # ---- Panel 2: rincian 3 sekuens terbaik (horizontal bar) ----
    label_top3 = [f"#{i+1}\n{rec['id'].split('|')[-1]}\n({rec['deskripsi'].split()[1]})"
                  for i, rec in enumerate(top3)]
    nilai_top3 = [rec["gc_content"] for rec in top3]
    y_pos = range(len(top3))

    bars2 = ax2.barh(y_pos, nilai_top3, color="#2E8B57", edgecolor="white")
    ax2.set_yticks(y_pos)
    ax2.set_yticklabels(label_top3, fontsize=9)
    ax2.invert_yaxis()
    ax2.set_xlabel("GC Content (%)", fontsize=10)
    ax2.set_title("3 Sekuens Terbaik\n(GC Content Tertinggi)", fontsize=12, fontweight="bold")
    ax2.set_xlim(0, max(nilai_top3) + 10)

    for bar, val in zip(bars2, nilai_top3):
        ax2.text(val + 0.8, bar.get_y() + bar.get_height() / 2, f"{val}%",
                  va="center", fontsize=10, fontweight="bold")

    plt.tight_layout()
    plt.savefig(output_path, dpi=150)
    plt.close()


def tulis_csv(data_terurut, output_path):
    """Menuliskan seluruh hasil analisis (terurut) ke file CSV."""
    with open(output_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([
            "Peringkat", "ID_GenBank", "Deskripsi", "Panjang_bp",
            "Jumlah_A", "Jumlah_T", "Jumlah_G", "Jumlah_C", "GC_Content_persen"
        ])
        for i, rec in enumerate(data_terurut, start=1):
            writer.writerow([
                i, rec["id"], rec["deskripsi"], rec["panjang"],
                rec["jumlah_A"], rec["jumlah_T"], rec["jumlah_G"], rec["jumlah_C"],
                rec["gc_content"]
            ])


if __name__ == "__main__":
    main()
