def processing(inputpath: str, outpath: str) -> str:
    with open(outpath, "w", encoding="utf-8") as fo:
        with open(inputpath, "r", encoding="utf-8") as f:
            for number, line in enumerate(f):
                line = line.strip()
                name, pseq, cseq, tai = line.split("\t")
                fo.write(f">{name}\n{cseq}\n")


if __name__ == "__main__":
    processing(
        "data/raw_data/hightTAI.tsv",
        "data/test_data/teach/high_all.fasta",
    )
