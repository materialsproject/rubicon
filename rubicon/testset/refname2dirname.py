import json

__author__ = 'xiaohuiqu'


from rubicon.testset.get_g3_benchmark import get_g3_bench_collection

def refname2dirname(mission):
    collection = get_g3_bench_collection()
    result_cursor = collection.find({"user_tags.mission": mission},
                             fields=['path', 'user_tags.fw_name'])
    calc_result = list(result_cursor)
    fw_name_2_path = {m['user_tags']['fw_name'].strip(): m['path'] for m in calc_result}

    with open("gauname2refname.json") as f:
        gau2ref = json.load(f)
    gau2ref.pop('NA')
    ref2gau = {v: k for k, v in gau2ref.items()}

    ref2path = {ref_name: fw_name_2_path[fw_name] for ref_name, fw_name in ref2gau.items() if fw_name in fw_name_2_path}

    return ref2path


if __name__ == '__main__':
    larry_ref2path = refname2dirname('G2-97 Test Set Benchmark (Larry Scheme)')
    with open('larry_ref2path.json', 'w') as f:
        json.dump(larry_ref2path, f, indent=4)
    shuyue_ref2path = refname2dirname('G2-97 Test Set Benchmark (Shyue Scheme)')
    with open('shuyue_ref2path.json', 'w') as f:
        json.dump(shuyue_ref2path, f, indent=4)