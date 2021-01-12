rule fetch_sample_list:
    input: config["input_dir"]
    output: out
    params:
        script=script_dir + "/generate_sample_list.py",
    script:
        "{params.script} {input}"