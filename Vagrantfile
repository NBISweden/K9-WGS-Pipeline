# -*- mode: ruby -*-
# vi: set ft=ruby :

work_image = '../nxf_work.vdi'

Vagrant.configure("2") do |config|
  config.vm.box = "ubuntu/xenial64"

  config.vm.provider "virtualbox" do |vb|
    vb.memory = "8128"
    vb.cpus = 4
    if ARGV[0] == "up" && !File.exist?(work_image)
        vb.customize ['createhd', '--filename', work_image,
                                  '--size', 50 * 1024,
                                  '--format', 'VDI',
                                  '--variant', 'fixed']
    end
    if ARGV[0] == "up"
        vb.customize ['storageattach', :id,
                        '--storagectl', "SCSI",
                        '--port', 2, '--device', 0,
                        '--type', 'hdd', '--medium', work_image]
    end
    if ARGV[0] == "destroy"
        vb.customize ['storageattach', :id,
                        '--storagectl', "SCSI",
                        '--port', 2, '--device', 0,
                        '--type', 'hdd', '--medium', 'none']
    end
  end

  config.vm.provision "shell", path: "scripts/vagrant-setup-disk.sh"
  config.vm.provision "shell", path: "scripts/vagrant-mount-disk.sh"
  config.vm.provision "shell", path: "scripts/vagrant-provision.sh"
end
