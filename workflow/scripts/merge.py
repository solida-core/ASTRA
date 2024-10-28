from subprocess import run
import sys

with open(snakemake.log[0], "w") as f:
    sys.stderr = sys.stdout = f
    print(len(snakemake.input))
    if len(snakemake.input) > 1:
        cmd = [snakemake.params['cmd'], 'merge',
               '--threads', str(snakemake.threads),
               '-O', snakemake.params['output_fmt']]
        if 'genome' in snakemake.params:
            cmd.append('--reference')
            cmd.append(snakemake.params['genome'])
        cmd.append("-o")
        out = "tmp.{}".format(snakemake.output[0])
        cmd.append(out)
        for i in snakemake.input:
            cmd.append(i)
        run(cmd)
        run(['mv', out, snakemake.output[0]])
    else:
        if snakemake.params['output_fmt'] == 'CRAM':
            _opt_ = '-C'
        else:
            _opt_ = '-b'
        cmd = [snakemake.params['cmd'], 'view', _opt_,
               '--threads', str(snakemake.threads),
               '-O', snakemake.params['output_fmt']]
        if snakemake.params['output_fmt'] == 'CRAM':
            cmd.append('--reference')
            cmd.append(snakemake.params['genome'])
        cmd.append('-o')
        out = "tmp.{}".format(snakemake.output[0])
        cmd.append(out)
        cmd.append(snakemake.input[0])
        run(cmd)
        run(['mv', out, snakemake.output[0]])
        run(['touch', '-h', snakemake.output[0]])
