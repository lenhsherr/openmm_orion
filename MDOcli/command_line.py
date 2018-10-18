import click

from openeye import oechem

from datarecord import read_mol_record

from Standards import Fields

from datarecord import (Types,
                        OEField,
                        OEWriteRecord)

from cuberecord import Types as TypesCR

import mdtraj as md

import os


@click.group(
    context_settings={
        "help_option_names": ("-h", "--help")
    }
)
@click.pass_context
def main(ctx):
    ctx.obj = dict()


@main.group()
@click.argument('filename', type=click.Path(exists=True))
@click.option("--id", help="Record ID number", default="all")
@click.option("--profile", help="OCLI profile name", default="default")
@click.pass_context
def dataset(ctx, filename, id, profile):
    """ Records Extraction"""

    ctx.obj['filename'] = filename

    ifs = oechem.oeifstream(filename)

    records = []

    while True:
        record = read_mol_record(ifs)
        if record is None:
            break
        records.append(record)
    ifs.close()

    if id == 'all':
        ctx.obj['records'] = records
    else:
        if int(id) < len(records):
            ctx.obj['records'] = [records[int(id)]]
        else:
            raise ValueError("Wrong record number selection: {} > max = {}".format(int(id), len(records)))

    os.environ['ORION_PROFILE'] = profile


@dataset.command("trajectory")
@click.option("--format", help="Trajectory format", type=click.Choice(['h5', 'dcd', 'tar.gz']), default='h5')
@click.option("--stgn", help="MD Stage number", default="last")
@click.option("--fixname", help="Edit the trajectory file name", default=None)
@click.pass_context
def trajectory_extraction(ctx, format, stgn, fixname):
    stagen = stgn

    new_record_list = []

    orion_trj_field = OEField("Trajectory_OPLMD", TypesCR.Orion.File)

    for record in ctx.obj['records']:

        if not record.has_value(Fields.md_stages):
            print("No MD stages have been found in the selected record")
            continue

        stages = record.get_value(Fields.md_stages)
        nstages = len(stages)
        title = record.get_value(Fields.title)
        id = record.get_value(Fields.id)
        fn = title + "_" + str(id)

        if stagen == 'last':
            stages_work_on = [stages[-1]]
        elif stagen == 'all':
            stages_work_on = stages
        else:
            if int(stagen) > nstages:
                print("Wrong stage number selection: {} > max = {}".format(stagen, nstages))
                return
            else:
                stages_work_on = [stages[int(stagen)]]

        for idx in range(0, len(stages_work_on)):

            stage = stages_work_on[idx]

            if stage.has_value(Fields.trajectory):
                fnt = stage.get_value(Fields.trajectory)

            elif stage.has_value(orion_trj_field):
                trj_id = stage.get_value(orion_trj_field)
                suffix = ''
                if stage.has_value(Fields.log_data):
                    log = stage.get_value(Fields.log_data)
                    log_split = log.split()
                    for i in range(0, len(log_split)):
                        if log_split[i] == 'suffix':
                            suffix = log_split[i+2]
                            break
                fnt = fn + '-' + suffix + '.' + format
                trj_id.download_to_file(fnt)
            else:
                print("No MD trajectory found in the selected stage record {}".format(stage.get_value(Fields.stage_name)))
                continue

            if format == 'dcd':
                print(fnt)
                trj_mdtraj = md.load(fnt)
                trj_mdtraj[0].save(os.path.splitext(fnt)[0]+'.pdb')
                trj_mdtraj.save(os.path.splitext(fnt)[0]+'.dcd')

            if fixname is not None:
                stage.set_value(Fields.trajectory, fnt)
                stages_work_on[idx] = stage

            record.set_value(Fields.md_stages, stages_work_on)
            new_record_list.append(record)

    if fixname is not None:

        ofs = oechem.oeofstream(fixname)

        for record in new_record_list:
            OEWriteRecord(ofs, record, fmt='binary')


@dataset.command("logs")
@click.option("--stgn", help="MD Stage number", default="last")
@click.pass_context
def logs_extraction(ctx, stgn):

    stagen = stgn

    for record in ctx.obj['records']:

        if not record.has_value(Fields.md_stages):
            print("No MD stages have been found in the selected record")
            return

        stages = record.get_value(Fields.md_stages)
        nstages = len(stages)
        title = record.get_value(Fields.title)
        id = record.get_value(Fields.id)
        fn = title + "_" + str(id)

        if stagen == 'last':
            stages_work_on = [stages[-1]]
        elif stagen == 'all':
            stages_work_on = stages
        else:
            if int(stagen) > nstages:
                print("Wrong stage number selection: {} > max = {}".format(stagen, nstages))
                return
            else:
                stages_work_on = [stages[int(stagen)]]

        for stage in stages_work_on:

            print("SYSTEM NAME = {}".format(fn))
            print("Stage name = {}".format(stage.get_value(Fields.stage_name)))
            if stage.has_value(Fields.log_data):
                print(stage.get_value(Fields.log_data))
            else:
                print("No logs have been found")


@dataset.command("protein")
@click.pass_context
def protein_extraction(ctx):

    for record in ctx.obj['records']:

        if not record.has_value(Fields.protein):
            print("No protein have been found in the selected record")
            return
        else:
            title = record.get_value(Fields.title).split("_")[0]
            id = record.get_value(Fields.id)
            fn = title + "_" + str(id)+".oeb"
            with oechem.oemolostream(fn) as ofs:
                oechem.OEWriteConstMolecule(ofs, record.get_value(Fields.protein))


@dataset.command("ligand")
@click.pass_context
def ligand_extraction(ctx):

    for record in ctx.obj['records']:

        if not record.has_value(Fields.ligand):
            print("No ligand have been found in the selected record")
            return
        else:
            title = record.get_value(Fields.title).split("_")[1:]
            title = "_".join(title)
            id = record.get_value(Fields.id)
            fn = title + "_" + str(id)+".oeb"
            with oechem.oemolostream(fn) as ofs:
                oechem.OEWriteConstMolecule(ofs, record.get_value(Fields.ligand))


@dataset.command("info")
@click.pass_context
def info_extraction(ctx):

    def recursive_record(record, level=0):

        for field in record.get_fields():
            if not issubclass(field.get_type(), Types.Record) or not issubclass(field.get_type(), Types.RecordVec):
                blank = (level+1) * "       "
                print("{} |".format(blank))
                print("{} |".format(blank))
                dis = "______"
                if (field.get_type() is Types.String or
                   field.get_type() is Types.Int or
                   field.get_type() is Types.Float):
                        print("{} {} name = {} type = {} value = {}".format(blank, dis,
                                                                            field.get_name(),
                                                                            field.get_type(),
                                                                            str(record.get_value(field))[0:30]))
                else:
                    print("{} {} name = {} type = {}".format(blank, dis,
                                                             field.get_name(),
                                                             field.get_type()))

            if issubclass(field.get_type(), Types.Record):
                recursive_record(record.get_value(field), level+1)
            elif issubclass(field.get_type(), Types.RecordVec):
                vec = record.get_value(field)
                for rec in vec:
                    recursive_record(rec, level+1)

    for record in ctx.obj['records']:
        recursive_record(record, 0)