from Bio import GenBank
import re as re
from pathlib import Path
from shiny import App, render, ui



def ui_card(title, *args):
    return (
        ui.div(
            {"class": "card mb-4"},
            ui.div(title, class_="card-header"),
            ui.div({"class": "card-body"}, *args),
        ),
    )
"\\wsl.localhost\Debian\home\sa\myapp"
app_ui = ui.page_fluid(
    ui.h2("Viral Segment Concatenator"),
    ui.input_file("file1", "Choose a file(.gb or fasta) to upload:", multiple=True),
    ui.input_switch("organism", "Name of your virus is contained in the field 'organism'"),
    ui.input_switch("only_proteint_coding", "Only protein-coding sequence"),
    ui.input_text("subtype", "Name of organism", "H5N5"),
    ui.input_numeric("nb_of_segments", "Number of segments", 8),
    ui.input_numeric("number_nn", "Number of characters (N) separating segments", 0),
    ui.input_numeric("nt", "Permissible difference", 200),
    ui.input_numeric("nb_segment1", "Segment 1", 2300),
    ui.input_numeric("nb_segment2", "Segment 2", 2300),
    ui.input_numeric("nb_segment3", "Segment 3", 2200),
    ui.input_numeric("nb_segment4", "Segment 4", 1700),
    ui.input_numeric("nb_segment5", "Segment 5", 1500),
    ui.input_numeric("nb_segment6", "Segment 6", 1400),
    ui.input_numeric("nb_segment7", "Segment 7", 1000),
    ui.input_numeric("nb_segment8", "Segment 8", 850),
    ui.output_text_verbatim("file_content"),
    ui_card(
        "The file is ready",
        ui.download_button("download1", "Download fasta"),
    ),
)


def server(input, output, session):

    @output
    
    @render.text
    def file_content():
        file_infos = input.file1()
        text = 'The number of sequences:'
        if not file_infos:
            return

        for file in file_infos:
            if file["datapath"].endswith('fasta'):
                subtype = input.subtype()
                organisms = read_fasta(file["datapath"], subtype)
            else:
                        
                organisms = {}
                
                #parsing of genbank
                with open(file["datapath"], "r") as handle:
                    
                    for record in GenBank.parse(handle):
                        

                        not_protein = []  
                        future_id_all = [""]
                        organism = ""
                        if input.organism():
                            organism = "/organism="
                        number_of_segment = 'wrong'
                        segments = [] 

                        for feature in record.features:
                                fut_id_var1 = ""
                                fut_id_var2 = ""
                                fut_id_var3 = ""
                                if feature.key == "CDS":
                                    not_protein.append(feature.location)
                                    
                                for qualifier in feature.qualifiers:  
                            
                                    if qualifier.key == "/segment=":  
                                        number_of_segment = qualifier.value.replace('"',"")
                                        segments.append(number_of_segment)
                                    if qualifier.key == "/strain=":  
                                        fut_id_var1 = qualifier.value
                                        future_id_all.append(fut_id_var1)
                                    if qualifier.key ==  "/isolate=":
                                        fut_id_var2 = qualifier.value
                                        future_id_all.append(fut_id_var2)
                                    #may be useful for some case
                                    if qualifier.key == organism:
                                        fut_id_var3 = qualifier.value
                                        future_id_all.append(fut_id_var3) 
                                    fut_id = max(future_id_all, key = len)
                                
                                if fut_id not in organisms.keys():
                                    organisms[fut_id] = {}
                        organisms[fut_id][number_of_segment] = [record.sequence]
                        if not_protein:
                            organisms[fut_id][number_of_segment].append(not_protein)
                            #organisms[fut_id][number_of_segment].append('..'.join(list(filter(None, re.split(r'\D',not_protein[0])))))   
                            #organisms[fut_id][number_of_segment].append('..'.join(list(filter(None, re.split(r'\D',not_protein[-1])))))    
                                        
            strings = []
            nb_of_segments = input.nb_of_segments()
            number_of_nn = input.number_nn()
            nt = input.nt()          
            nb_segments = {
                '1':  input.nb_segment1(),
                '2':  input.nb_segment2(),
                '3':  input.nb_segment3(),
                '4':  input.nb_segment4(),
                '5':  input.nb_segment5(),
                '6':  input.nb_segment6(),
                '7':  input.nb_segment7(),
                '8':  input.nb_segment8(),
            }
            
            organisms = change_letters_to_numbers(organisms)
            organisms = choose_the_proteint_coding_sequence(organisms)
            if input.only_proteint_coding():
                organisms = cut_primers_for_all_segments(organisms)
            else:
                organisms = cut_primers(organisms, nb_of_segments)
            organisms_right_len, t = check_segments(organisms, nb_segments, nt, number_of_nn)
            prep_without_dn = prep_delete_degenerate_nucleotides(organisms_right_len, nb_of_segments)
            without_dn = delete_degenerate_nucleotides(prep_without_dn)

            #create a file with mergening genome of virus
            with open("all_segments.fasta", "w") as f:
                for key1, item in  without_dn.items():
                    strings.append(">{} \n{}\n".format(key1, item))
                result1 = "".join(strings)
                f.write(result1)
            
       
        return text, len(without_dn)
    @session.download()
    def download1(): 
        path = Path(__file__).parent / "all_segments.fasta"
        
        return str(path)
    
#read fasta
def read_fasta(file, subtype):
    segments = {
    'PB2': "1",
    'hemagglutinin': '4',
    'neuraminidase': '6',
    
    "PB2": "2",
    "PA": "3",
    "HA": "4",
    "NP": "5",
    "NA": "6",
    "M1": "7",
    "NS1": "8"
        
    }
    with open(file) as fasta:
        organisms = {}
        name =""
        v_seg = 0
        flag = False
        for line in fasta:
            if line.startswith('>'):
                name = line.split(subtype)[0]
                if len(name) >= 2 and '(' in name:
                    name = name.split('(')[1]
                
            
            
                if "segment" in line:
                    v_seg = line.split("segment")[1][1]
                    if name not in organisms.keys():
                        organisms[name] = {}
                    organisms[name][v_seg] = ["",[""]]
                    flag = True
                elif "gene" in line:
                    v_seg = line.split("gene")[0].split()[-1].strip("()")
                    if v_seg in segments.keys():
                        v_seg = segments[v_seg]
                        if name not in organisms.keys():
                            organisms[name] = {}
                        organisms[name][v_seg] = ["",[""]]
                        flag = True
                    else:
                        flag = False
                else:
                    flag = False
            elif flag:
                organisms[name][v_seg][0]+=line.strip()
                organisms[name][v_seg][1][0] = f'1..{len(organisms[name][v_seg][0])}'
            
    return organisms


#in some cases, segment names consist of letters, we need to change them to numbers
def change_letters_to_numbers(organisms):
    for seg_seq in organisms.values():
            not_allowed = []
            if 'M' in seg_seq.keys() and 'L' in seg_seq.keys() and 'S' in seg_seq.keys():
                    seg_seq['1']=seg_seq['S']
                    seg_seq.pop('S')
                    seg_seq['2']=seg_seq['M']
                    seg_seq.pop('M')
                    seg_seq['3']=seg_seq['L']
                    seg_seq.pop('L')
            elif 'S' in seg_seq.keys() and 'L' in seg_seq.keys() and 'M' not in seg_seq.keys():
                    seg_seq['1']=seg_seq['S']
                    seg_seq.pop('S')
                    seg_seq['2']=seg_seq['L']
                    seg_seq.pop('L')
            not_allowed = list(filter(None, "".join(re.split(r'\d', "/".join(seg_seq.keys()))).replace('"',"").split("/")))
            
            if len(not_allowed) != 0:
                for wrong_name in not_allowed:
                      seg_seq.pop(wrong_name)
           
    return organisms
#we have to choose the longest segment

def choose_the_proteint_coding_sequence(lib):
    for name, content in lib.items():
        for num, data in content.items():
            #print(data)
            if len(data) > 1:
                
                
                good = []
                for item in data[1]:
                    if "(" in item:
                        good.extend(item.split("(")[1][:-1].split(","))
                    else:
                        good.append(item)
                max = 0
                best = ""
                #print(good)
                for item in good:
                    t = item.split("..")
                    #print(item)
                    i = int(t[0].strip("<>"))
                    j = int(t[1].strip("<>"))
                    if j - i > max:
                        max = j-i
                        best = item
                
                data[1] = best.replace(">","").replace("<","")
                
            else:
                data = []
    
    return lib
                

#sometimes we should cut nucleotides before first start-codon and after the last stop-codon
def cut_primers(organisms,  nb_of_segments):

    for dictt in organisms.values():
        for number_of_segment, listt in dictt.items():

            if number_of_segment == '1':
                if len(listt) > 1:
                    cut_start = listt[-1].split("..")[0]
                    #if cut_start not in ["A", "T", "G", "C", "N"]:
                    cut_start = int(cut_start) - 1

                    dictt[number_of_segment] = listt[0][cut_start::]
               
            if int(number_of_segment.replace('"', "")) == nb_of_segments:
                if len(listt) > 1:
                    cut_end = listt[-1].split("..")[-1]
                    #if cut_end not in ["A", "T", "G", "C"]:
                    cut_end = int(cut_end)
                    dictt[number_of_segment] = listt[0][0:cut_end]

            elif number_of_segment !='1' and number_of_segment != nb_of_segments:
                  dictt[number_of_segment] = listt[0]
                     
    return organisms



def cut_primers_for_all_segments(organisms):

    for dictt in organisms.values():
        for number_of_segment, listt in dictt.items():

                if len(listt) > 1:
                    cut_start = listt[-1].split("..")[0]
                    cut_start = int(cut_start) - 1
                    cut_end = listt[-1].split("..")[-1]
                    cut_end = int(cut_end)
                    
                    dictt[number_of_segment] = listt[0][cut_start:cut_end]

                     
    return organisms


#check len of segments
def check_segments(organisms, nb_segments, nt, number_of_nn):
    seq2 = {}
    t = {}
    for name, dictt in organisms.items():
        seq2[name] = {}  
        t[name] = {}
        for k, v in nb_segments.items():  
             if k in dictt.keys():
                t[name][k] = len(dictt[k]) 
                if (v - nt <= len(dictt[k]) <= v + nt):
                             seq2[name][k]=''.join(dictt[k]) + number_of_nn*"-"
                             

    return seq2, t                 

#here we bring to the form {id:record of all segments}
def prep_delete_degenerate_nucleotides(seq2, nb_of_segments):
    prep_without_dn = {}
    for key, values in seq2.items():
        if len(values) == nb_of_segments:
         prep_without_dn[key] = "".join(str(dict(sorted(values.items(), key=lambda x: x[0])).values())).replace("'","").replace('dict_values([',"").replace('])',"").replace(',',"").replace(" ","")
         #prep_without_dn[key] = dict(sorted(values.items(), key=lambda x: x[0])).values()
   
    return prep_without_dn

#sometimes we should delete {id:record} with degenerate nucleotides
def delete_degenerate_nucleotides(prep_without_dn):
    without_dn={}
    dn = ['u', 'w', 'n', 's', 'm', 'k', 'r', 'y', 'b', 'd', 'h', 'v']
    for key3, value3 in prep_without_dn.items():
        y = 0
        for nt in dn:
            if nt in value3.lower():
                y = 1
        if y == 0:
            without_dn[key3]=value3
    return without_dn





app = App(app_ui, server)