from subprocess import run
import sys

with open(snakemake.log[0], "w") as f:
    sys.stderr = sys.stdout = f
    print(len(snakemake.input['bams']))
    if len(snakemake.input['bams']) > 1:
        cmd = ['samtools', 'merge', '--threads', str(snakemake.threads), '-O', snakemake.params['output_fmt'],
               snakemake.output['cram']]
        # cmd.append("-o")
        for i in snakemake.input['bams']:
            cmd.append(i)
        run(cmd)

        cmd = ['samtools', 'index',
               '--threads', str(snakemake.threads),
               snakemake.output['cram']]
        run(cmd)
    else:
        if snakemake.params['output_fmt'] == 'CRAM':
            _opt_ = '-C'
        else:
            _opt_ = '-b'
        cmd = [snakemake.params['cmd'], 'view', _opt_, '--threads', str(snakemake.threads), '-O',
               snakemake.params['output_fmt'], snakemake.output['cram'], snakemake.input['bams'].pop()]

        # cmd.append('-o')
        run(cmd)

        run(['touch', '-h', snakemake.output[0]])

        cmd = ['samtools', 'index',
               '--threads', str(snakemake.threads),
               snakemake.output['cram']]
        run(cmd)
