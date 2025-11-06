"""
import_triage_dag.py
Imports Triage study to MySQL database.
"""
import os
import sys

from airflow.models.param import Param
from airflow.decorators import task
from airflow.exceptions import AirflowException
from airflow.models import DagRun
from airflow.utils.state import DagRunState

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from dags.import_base import ImporterConfig, build_import_dag

@task
def check_tempo_dag_not_running(**context):
    """Check if import_tempo_dag is running. If so, fail this DAG run."""
    from airflow.models import DagRun
    from airflow.utils.state import DagRunState

    active_runs = DagRun.find(
        dag_id="import_tempo_dag",
        state=DagRunState.RUNNING,
    )

    if active_runs:
        raise AirflowException(
            f"import_tempo_dag is currently running. "
            f"Cannot start import_triage_dag until import_tempo_dag completes."
        )

    return True

def _wire(tasks: dict[str, object]) -> None:
    check_task = check_tempo_dag_not_running()
    check_task >> tasks["data_repos"]
    tasks["data_repos"] >> tasks["fetch_data"]
    tasks["fetch_data"] >> tasks["setup_import"]
    tasks["setup_import"] >> tasks["import_sql"] >> tasks["cleanup_data"]

_TRIAGE_CONFIG = ImporterConfig(
    dag_id="import_triage_dag",
    description="Imports Triage study to MySQL database",
    importer="triage",
    tags=["triage"],
    target_nodes=("pipelines3_ssh",),
    data_nodes=("pipelines3_ssh",),
    task_names=(
        "fetch_data",
        "setup_import",
        "import_sql",
        "cleanup_data",
    ),
    db_properties_filename="manage_triage_database_update_tools.properties",
    color_swap_config_filename=None, # Not used for MySQL
    params={
        "data_repos": Param(
            [
                "datahub",
                "cmo-argos",
                "private",
                "impact"
            ],
            type="array",
            description="Comma-separated list of data repositories to pull updates from/cleanup.",
            title="Data Repositories",
            examples=[
                "datahub",
                "bic-mskcc-legacy",
                "cmo-argos",
                "private",
                "impact",
                "knowledge-systems-curated-studies",
                "datahub_shahlab",
                "msk-mind-datahub",
                "pipelines-testing",
                "genie",
                "extract-projects",
                "cmo-access"
            ],
        ),
    },
    wire_dependencies=_wire,
    pool="triage_import_pool",
    schedule_interval="0 0 * * *",
)

globals()[_TRIAGE_CONFIG.dag_id] = build_import_dag(_TRIAGE_CONFIG)
