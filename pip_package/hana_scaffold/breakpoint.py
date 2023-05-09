# -*- coding: utf-8 -*-
import json
import subprocess
import os.path
import sys

_BP_INIT = False
_BP_PATH = ''


class BreakPoint(dict):
    def __init__(self):
        super().__init__()
        self['finished'] = False
        self['step_last'] = -1
        self['step_current'] = -1

    def load(self):
        try:
            with open(_BP_PATH, 'r') as bp_f:
                self.update(json.load(bp_f))
        except Exception:
            raise Exception('Bad breakpoint detected.')
        # Reset the step current to first.
        if _BP_DATA['finished'] or '-f' in sys.argv:
            _BP_DATA['finished'] = False
            self['step_last'] = -1
        else:
            self['step_last'] = self['step_current']
        self['step_current'] = -1

    def save(self):
        global _BP_INIT, _BP_PATH
        if _BP_INIT:
            with open(_BP_PATH, 'w') as bp_f:
                json.dump(self, bp_f, indent='\t')


_BP_DATA = BreakPoint()


def init():
    global _BP_INIT, _BP_PATH, _BP_DATA
    if _BP_INIT:
        raise Exception('Breakpoint is already initialized.')
    script_path = sys.argv[0]
    script_dir, _ = os.path.split(script_path)
    _BP_PATH = os.path.join(script_dir, '.hana_breakpoint')
    if os.path.isfile(_BP_PATH):
        # Load the breakpoint info.
        _BP_DATA.load()
    _BP_INIT = True


def finish():
    # Save the breakpoint out.
    global _BP_INIT, _BP_DATA
    _BP_DATA['finished'] = True
    if _BP_INIT:
        _BP_DATA.save()


def _breakpoint_execute(callback):
    global _BP_DATA
    # Update the current step ID.
    _BP_DATA['step_current'] += 1
    # If breakpoint is initialized and passed the last step.
    if not _BP_INIT or _BP_DATA['step_current'] > _BP_DATA['step_last']:
        # This step cannot be skipped.
        callback()
    # Save the breakpoint.
    _BP_DATA.save()


def run_command(mod_args: list):
    def execute():
        print('> {}'.format(' '.join(mod_args)))
        mod_proc = subprocess.Popen(mod_args)
        mod_proc.wait()
    _breakpoint_execute(execute)


def run_pipeline(mod_pipeline: list):
    if len(mod_pipeline) == 1:
        return run_command(mod_pipeline[0])

    def execute():
        print('>', ' | '.join(' '.join(x) for x in mod_pipeline))
        # Run the first command.
        last_proc = subprocess.Popen(mod_pipeline[0], stdout=subprocess.PIPE)
        # Run the body commands.
        for mod_body in mod_pipeline[1:-2]:
            last_proc = subprocess.Popen(mod_body, stdin=last_proc.stdout, stdout=subprocess.PIPE)
        # Run the last command.
        last_proc = subprocess.Popen(mod_pipeline[-1], stdin=last_proc.stdout)
        last_proc.wait()
    _breakpoint_execute(execute())
