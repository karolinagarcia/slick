parent_path=$( cd "$(dirname "${BASH_SOURCE[0]}")" ; pwd -P )
chmod u+x "$parent_path/bin/slick_frontend.sh"
chmod u+x "$parent_path/bin/slick_new.sh"
echo "export PATH=\"$parent_path/bin:\$PATH\"" >> $HOME/.bashrc
