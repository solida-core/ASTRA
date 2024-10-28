from subprocess import run
import sys

with open(snakemake.log[0], "w") as f:
    sys.stderr = sys.stdout = f
    print(len(snakemake.input['bams']))
    if len(snakemake.input['bams']) > 1:
        cmd = ['samtools', 'merge',
               '--threads', str(snakemake.threads),
               '-O', snakemake.params['output_fmt']]
        if 'genome' in snakemake.params:
            cmd.append('--reference')
            cmd.append(snakemake.params['genome'])
        cmd.append("-o")
        cmd.append(snakemake.output['cram'])
        for i in snakemake.input['bams']:
            cmd.append(i)
        run(cmd)

        cmd = ['samtools', 'index',
               '--threads', str(snakemake.threads),
               '--output', snakemake.output['crai'],
               snakemake.output['cram']]
        run(cmd)
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
        cmd.append(snakemake.output['cram'])
        cmd.append(snakemake.input['bams'].pop())
        run(cmd)

        run(['touch', '-h', snakemake.output[0]])

        cmd = ['samtools', 'index',
               '--threads', str(snakemake.threads),
               '--output', snakemake.output['crai'],
               snakemake.output['cram']]
        run(cmd)
