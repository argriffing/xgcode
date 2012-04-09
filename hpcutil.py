"""
This module is for running remote jobs.

It is designed for high-cpu and small-data batches,
with no custom code on the server side.
"""

from StringIO import StringIO
import subprocess
import os
import time

import Util

def makedirs_checked(s):
    if not os.path.exists(s):
        os.makedirs(s)

class Location:
    def __init__(self, batch_root, name):
        self.batch_root = batch_root
        self.name = name
    def get_batch_in(self):
        return os.path.join(self.batch_root, 'batch_in')
    def get_batch_out(self):
        return os.path.join(self.batch_root, 'batch_out')
    def get_in(self):
        return os.path.join(self.get_batch_in(), 'in_' + self.name)
    def get_out(self):
        return os.path.join(self.get_batch_out(), 'out_' + self.name)
    def get_in_bsubs(self):
        return os.path.join(self.get_in(), 'bsubs')
    def get_in_contents(self):
        return os.path.join(self.get_in(), 'contents')
    def get_in_tgz(self):
        return self.get_in() + '.tgz'
    def get_out_tgz(self):
        return self.get_out() + '.tgz'

class RemoteBase:
    """
    Supply the preprocess member function.
    The getp prefixed functions get paths to files and to directories.
    """
    def __init__(
            self,
            local_batch_root,
            remote_batch_root,
            remote_server):
        self.name = self._reserve_name()
        self.remote_server = remote_server
        self.local = Location(local_batch_root, self.name)
        self.remote = Location(remote_batch_root, self.name)
    def _remote_lsf_call(self, args):
        full_args = [
                'source',
                '/usr/local/lsf/conf/cshrc.lsf',
                ';'] + list(args)
        return self._remote_call(tuple(full_args))
    def _remote_call(self, args):
        full_args = ['ssh', self.remote_server] + list(args)
        return subprocess.check_output(
                tuple(full_args),
                stderr=subprocess.STDOUT)
    def _reserve_name(self):
        out = StringIO()
        print >> out, 'reserving a unique random local name...'
        data = out.getvalue()
        tmp_path = Util.create_tmp_file(data=data, prefix='', suffix='')
        tmp_name = os.path.basename(tmp_path)
        t = time.gmtime()
        prefix = '%04d%02d%02d_' % (t.tm_year, t.tm_mon, t.tm_mday)
        return prefix + tmp_name
    def _init_local_structure(self):
        makedirs_checked(self.local.get_in_bsubs())
        makedirs_checked(self.local.get_in_contents())
        makedirs_checked(self.local.get_batch_out())
    def _make_local_tgz(self):
        tgz_local = self.local.get_in_tgz()
        args = ('tar', 'czvf', tgz_local,
                '-C', self.local.get_batch_in(), 'in_' + self.name)
        result = subprocess.check_output(args, stderr=subprocess.STDOUT)
    def _init_remote_structure(self):
        for d in (self.remote.get_batch_in(), self.remote.get_out()):
            args = ('mkdir', '-p', d)
            self._remote_call(args)
    def _scp_to_remote(self):
        tgz_local = self.local.get_in_tgz()
        tgz_remote = self.remote.get_in_tgz()
        args = ('scp', tgz_local, '%s:%s' % (self.remote_server, tgz_remote))
        result = subprocess.check_output(args, stderr=subprocess.STDOUT)
    def _remote_expand(self):
        tgz_remote = self.remote.get_in_tgz()
        args = ('tar', 'xzvf', tgz_remote, '-C', self.remote.get_batch_in())
        result = self._remote_call(args)
    def _submit(self):
        for name in os.listdir(self.local.get_in_bsubs()):
            if name.endswith('.bsub'):
                remote_path = os.path.join(self.remote.get_in_bsubs(), name)
                args = ('bsub', '<', remote_path)
                result = self._remote_lsf_call(args)
    def _poll(self):
        while True:
            args = ('bjobs',)
            result = self._remote_lsf_call(args)
            if result:
                lines = [line.strip() for line in result.splitlines()]
                lines = [line for line in lines if line]
                print result
                if lines[-1] == 'No unfinished job found':
                    print 'umm no jobs are left'
                    break
            print 'jobs remain (?)'
            time.sleep(2)
    def _make_remote_tgz(self):
        args = ('tar', 'czvf', self.remote.get_out_tgz(),
                '-C', self.remote.get_batch_out(), 'out_' + self.name)
        result = self._remote_call(args)
    def _scp_to_local(self):
        args = ('scp',
                '%s:%s' % (self.remote_server, self.remote.get_out_tgz()),
                self.local.get_out_tgz())
        result = subprocess.check_output(args, stderr=subprocess.STDOUT)
    def _local_expand(self):
        args = ('tar', 'xzvf', self.local.get_out_tgz(),
                '-C', self.local.get_batch_out())
        result = subprocess.check_output(args, stderr=subprocess.STDOUT)
    def run(self):
        """
        This is a long blocking call that runs the batch remotely.
        It expects that ssh keys have been set up.
        """
        # (local) reserve a unique identifier
        self._reserve_name()
        # (local) set up the local directory structure
        self._init_local_structure()
        # (local) call preprocess callback to populate batch_in
        self.preprocess()
        # (local) tgz batch_in
        self._make_local_tgz()
        # (ssh) set up the remote directory structure
        self._init_remote_structure()
        # (scp) scp the tgz to the remote server
        self._scp_to_remote()
        # (ssh) expand the tgz
        self._remote_expand()
        # (ssh, lsf) submit all of the bsub files
        self._submit()
        # (ssh, lsf) poll for completion of the jobs
        self._poll()
        # (ssh) tgz the results
        self._make_remote_tgz()
        # (scp) scp the tgz to the local computer
        self._scp_to_local()
        # (local) expand the tgz
        self._local_expand()


class RemoteBrc(RemoteBase):
    def __init__(self):
        RemoteBase.__init__(
                self,
                '/tmp',
                '/brc_share/brc/argriffi',
                'login02.hpc.ncsu.edu')

class RemoteBrcSilly(RemoteBrc):
    def preprocess(self):
        local_dummy_path = os.path.join(
                self.local.get_in_contents(), 'hello.world')
        remote_dummy_path = os.path.join(
                self.remote.get_in_contents(), 'hello.world')
        with open(local_dummy_path, 'w') as fout:
            print >> fout, 'hello world'
            print >> fout, 'goodbye world'
        # define the first bsub job
        bsub_path = os.path.join(self.local.get_in_bsubs(), 'job1.bsub')
        remote_result_path = os.path.join(self.remote.get_out(), 'foo.txt')
        with open(bsub_path, 'w') as fout:
            print >> fout, 'wc %s > %s' % (
                    remote_dummy_path, remote_result_path)
        # define the second bsub job
        bsub_path = os.path.join(self.local.get_in_bsubs(), 'job2.bsub')
        remote_result_path = os.path.join(self.remote.get_out(), 'bar.txt')
        with open(bsub_path, 'w') as fout:
            print >> fout, 'wc -l %s > %s' % (
                    remote_dummy_path, remote_result_path)


if __name__ == '__main__':
    RemoteBrcSilly().run()

