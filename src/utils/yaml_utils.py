import yaml


def load_yaml(yaml_fn):
    with open(yaml_fn, encoding="UTF-8") as handler:
        return yaml.safe_load(handler)


def save_yaml(yaml_fn, data):
    with open(yaml_fn, "w", encoding="UTF-8") as handler:
        yaml.dump(data, handler, default_flow_style=False)
