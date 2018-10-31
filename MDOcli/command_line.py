import click

from openeye import oechem

from datarecord import read_mol_record

from Standards import Fields

from datarecord import (Types,
                        OEWriteRecord)

from datarecord.datarecord import RecordVecData, RecordData

from orionclient.types import File

import mdtraj as md

import os

from orionclient.session import (OrionSession,
                                 get_profile_config,
                                 get_session)


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
def dataset(ctx, filename, id, profile=None, max_retries=2):
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
    # TODO
    if profile == "default" and os.environ.get("ORION_PROFILE") is not None:
        profile = os.environ["ORION_PROFILE"]

    profile_config = get_profile_config(profile=profile)

    ctx.obj['profile'] = profile
    ctx.obj['session'] = OrionSession(
            config=profile_config,
            requests_session=get_session({404: max_retries}))

    ctx.obj['credentials'] = profile_config

    if ctx.obj["credentials"] is not None:
        ctx.obj['session'].config = ctx.obj["credentials"]
    else:
        click.secho("Unable to find credentials", fg='red', err=True)


@dataset.command("trajectory")
@click.option("--format", help="Trajectory format", type=click.Choice(['h5', 'dcd', 'tar.gz']), default='h5')
@click.option("--stgn", help="MD Stage number", default="last")
@click.option("--fixname", help="Edit the trajectory file name", default=None)
@click.pass_context
def trajectory_extraction(ctx, format, stgn, fixname):
    stagen = stgn

    new_record_list = []

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

            elif stage.has_value(Fields.orion_local_trj_field):

                trj_id = stage.get_value(Fields.orion_local_trj_field)

                suffix = ''
                if stage.has_value(Fields.log_data):
                    log = stage.get_value(Fields.log_data)
                    log_split = log.split()
                    for i in range(0, len(log_split)):
                        if log_split[i] == 'suffix':
                            suffix = log_split[i+2]
                            break
                fnt = fn + '-' + suffix + '.' + format

                resource = ctx.obj['session'].get_resource(File, trj_id)

                resource.download_to_file(fnt)

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
            print("Stage type = {}".format(stage.get_value(Fields.stage_type)))
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

            field_type = field.get_type()

            blank = "       "
            print("{} |".format(blank * (level + 1)))
            print("{} |".format(blank * (level + 1)))
            dis = "______"

            if not field_type == RecordData and not field_type == RecordVecData:

                if (field.get_type() is Types.String or
                        field.get_type() is Types.Int or
                        field.get_type() is Types.Float):
                    print("{} {} name = {} type = {} value = {}".format(blank * (level + 1), dis,
                                                                        field.get_name(),
                                                                        field.get_type(),
                                                                        str(record.get_value(field))[0:15]))
                else:
                    print("{} {} name = {} type = {}".format(blank * (level + 1), dis,
                                                             field.get_name(),
                                                             field.get_type()))

            elif field_type == RecordData:
                print("{} {} RECORD".format(blank * (level + 1), dis))
                recursive_record(record.get_value(field), level + 1)

            elif field_type == RecordVecData:
                vec = record.get_value(field)
                print("{} {} RECORD VECTOR containing {} records".format(blank * (level + 1), dis, len(vec)))
                print("{} |".format(blank * (level + 2)))
                print("{} |".format(blank * (level + 2)))

                for idx in range(0, len(vec)):
                    print("{} {} RECORD # {}".format(blank * (level + 2), dis, idx))
                    recursive_record(vec[idx], level + 2)
                    if idx != len(vec) - 1:
                        print("{} |".format(blank * (level + 2)))
                        print("{} |".format(blank * (level + 2)))
            else:
                raise ValueError("Field type error: {}".format(field_type))

    for idx in range(0, len(ctx.obj['records'])):
        print(30 * "*" + " RECORD {}/{} ".format(idx + 1, len(ctx.obj['records'])) + 30 * "*")
        recursive_record(ctx.obj['records'][idx], 0)
        print("\n" + 30 * "*" + " END RECORD ".format(idx + 1, len(ctx.obj['records'])) + 30 * "*" + "\n") 