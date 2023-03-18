# -*- coding: utf-8 -*-
import json
import os.path

TEMPLATES = {
    'configure': '{}.conf.json',
    'status': '{}.status.json'
}


def project_path(dir_path: str, template_id: str, project_id: str):
    return os.path.join(dir_path, TEMPLATES[template_id].format(project_id))


def load_json(dir_path: str, template_id: str, project_id: str):
    try:
        with open(project_path(dir_path, template_id, project_id), 'r', encoding='utf-8') as json_file:
            return json.load(json_file)
    except Exception:
        return {}


def save_json(content, dir_path: str, template_id: str, project_id: str):
    with open(project_path(dir_path, template_id, project_id), 'w', encoding='utf-8') as json_file:
        return json.dump(content, json_file)


def load_project(dir_path: str, project_id: str):
    return load_json(dir_path, 'configure', project_id)


def load_progress(dir_path: str, project_id: str):
    return load_json(dir_path, 'status', project_id)
