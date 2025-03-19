# for generate with rarecodon

from accelerate import Accelerator
from Bio import SeqIO
from Bio.Seq import Seq

from DeepCodon.config.config import *
from DeepCodon.config.server_config import *
from DeepCodon.src.generate.beam import beam_search
from DeepCodon.src.generate.convert import index2content
from DeepCodon.src.generate.greedy import greedy_algorithm
from DeepCodon.src.generate.mixture import mixture_algorithm
from DeepCodon.src.generate.temperature import temperature_algorithm
from DeepCodon.src.generate.top_k import top_k_algorithm
from DeepCodon.src.generate.top_p import top_p_algorithm
from DeepCodon.src.services.preGpuQuery import *
from DeepCodon.src.tests.tai.centerCont import *
from DeepCodon.src.train.datasets import *
from DeepCodon.src.train.model import TransformerTranslator


def modelProcess(
    nowUuid: str,
    fastafile: str,
    smoothTokenDict: dict,
    generatenum: int,
    output_path: str,
    logging: object,
) -> str:
    if output_path is not None:
        outpath = output_path
    else:
        fastafilename = fastafile.split("/")[-1]
        print(f"fastafilename: {fastafilename}")
        outpath = f"{defeat_output_path}/{fastafilename}_{nowUuid}_{ifcase}_{generationmethod}.csv"
        outcachepath = f"{model_cache_path}/{fastafilename}_{nowUuid}_{ifcase}_{generationmethod}.csv"

    selected_gpu = assign_gpu(2048)
    if selected_gpu is not None:
        print(f"to GPU: {selected_gpu}")
    else:
        print("no GPU")
    logging.info(f"to GPU: {selected_gpu}")
    torch.cuda.set_device(selected_gpu)

    accelerator = Accelerator()
    device = accelerator.device
    with open(f"{model_cache_path}/{nowUuid}.tsv", "w") as tsv_file:
        for record in tqdm(SeqIO.parse(fastafile, "fasta")):
            # dna or protein
            inputseq = str(record.seq)
            if "M" not in inputseq:
                # Determine whether it can be divided by 3
                whichdata = "codon"
                if len(inputseq) % 3 != 0:
                    continue
                else:
                    protein_sequence = str(Seq(inputseq).translate())
            else:
                whichdata = "protein"
                protein_sequence = inputseq
            if len(protein_sequence) <= pclen and len(protein_sequence) >= minlen:
                if generatenum > 1:
                    for i in range(generatenum):
                        tsv_file.write(
                            f"{record.id}_{i}\t{protein_sequence}\t{inputseq}\n"
                        )
                else:
                    tsv_file.write(f"{record.id}\t{protein_sequence}\t{inputseq}\n")
    record_ids, enc_inputs, inputseqs, opmask_inputs = make_data_user(
        f"{model_cache_path}/{nowUuid}.tsv"
    )
    user_data = MyDataSetuser(record_ids, enc_inputs, inputseqs, opmask_inputs)
    user_loader = Data.DataLoader(user_data, shuffle=False)
    model = TransformerTranslator()
    model, eval_dataloader = accelerator.prepare(model, user_loader)
    model = accelerator.unwrap_model(model)
    model.load_state_dict(torch.load(whichmodel2generate, map_location=device))
    translations = []
    for number, (record_ids, enc_inputs, inputseqs, opmask_inputs) in tqdm(
        enumerate(user_loader)
    ):
        print(f"record_ids: {record_ids}\n")
        if number < maxgenerateseq:
            # Output model results
            enc_inputs = enc_inputs.cuda()
            opmask_inputs = opmask_inputs.cuda()
            if use_logit_bias:
                smooth_token_info = smoothTokenDict.get(record_ids[0])
                if smooth_token_info:
                    print("found smoothTokenDict")
                    logit_bias = smooth_token_info[3]
                else:
                    logit_bias = None
            else:
                logit_bias = None
            if generationmethod == "greedy":
                predict = greedy_algorithm(
                    model=model,
                    device=device,
                    src=enc_inputs,
                    opmask_inputs=opmask_inputs,
                )
            elif generationmethod == "beam":
                predict = beam_search(
                    model=model,
                    device=device,
                    src=enc_inputs,
                    opmask_inputs=opmask_inputs,
                )
            elif generationmethod == "top_p":
                predict = top_p_algorithm(model=model, device=device, src=enc_inputs)
            elif generationmethod == "top_k":
                predict = top_k_algorithm(model=model, device=device, src=enc_inputs)
            elif generationmethod == "temperature":
                predict = temperature_algorithm(
                    model=model, device=device, src=enc_inputs
                )
            elif generationmethod == "mix":
                # TODO
                predict = mixture_algorithm(
                    model=model,
                    device=device,
                    src=enc_inputs,
                    opmask_inputs=opmask_inputs,
                    debugFileOut=f"data/userData/output/debug/{number}",
                    logit_bias=logit_bias,
                )
            else:
                break
            output_str = index2content(
                enc_inputs, predict, inputseqs, record_ids, whichdata
            )
            translations.append(output_str)
        else:
            break

    print("finish genereate")
    with open(
        outpath,
        "w",
        encoding="utf-8",
    ) as f:
        if whichdata == "codon":
            f.write(
                "id,inputseq,outputseq,inputMFE,outputMFE,inputGC,outputGC,inputCAI,outputCAI,istrue\n"
            )
        else:
            f.write("id,inputseq,outputseq,outputMFE,outputGC,outputCAI,istrue\n")
        for translation in translations:
            f.write(translation + "\n")

    if whichdata == "codon":
        process = ProcessData(
            filepath=outpath,
            inputCol=1,
            outputCol=2,
        )
        process.processing(outpath=outcachepath)
        dfnotai = pd.read_csv(outpath)
        dftai = pd.read_csv(outcachepath)
        dftai = dftai.T
        dftai.reset_index(drop=True, inplace=True)
        dftai.drop(index=0, inplace=True)
        dftai.reset_index(drop=True, inplace=True)
        dftai.columns = ["inputTAI", "outputTAI"]
        df = pd.concat([dfnotai, dftai], axis=1)
        df.to_csv(
            outpath,
            index=False,
        )
    return outpath
