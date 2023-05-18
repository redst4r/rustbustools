use std::collections::VecDeque;
// use std::sync::mpsc::{channel, Sender, Receiver, sync_channel, SyncSender};
use crossbeam;
use crossbeam::channel::{unbounded, bounded, SendError, RecvError};

use std::thread;



pub struct BufferedSender<T>{
    buffer: VecDeque<T>,
    buffersize: usize,
    tx: crossbeam::channel::Sender<VecDeque<T>> // the inner sender has to handle collections of T!!
}

impl<T> BufferedSender<T>{
    pub fn new(tx: crossbeam::channel::Sender<VecDeque<T>>, buffersize: usize)-> Self{
        let buffer = VecDeque::with_capacity(buffersize); // note: the buffer has variable size!can get larger than buffersize
        BufferedSender {buffer,buffersize, tx}
    }

    pub fn send(&mut self, item: T) -> Result<(), SendError<VecDeque<T>>>{
        // sending a single item across the buffer
        // note that its not actually sent immediately, only once the buffer is full
        // or a .flush is triggered
        self.buffer.push_back(item);

        // if the buffer is full, send its content, clean the buffer to empty
        if self.buffer.len() >= self.buffersize{
            self.flush()
        }
        else{
            Ok(())
        }
    }
    fn flush(&mut self) -> Result<(), SendError<VecDeque<T>>>{
        println!("Flushing buffer -> sending all");
        let newbuffer: VecDeque<T> = VecDeque::with_capacity(self.buffersize);
        // swapping out the current buffer for an empty one
        let current_buf = std::mem::replace(&mut self.buffer,  newbuffer);
        self.tx.send(current_buf)
    }
}

pub struct BufferedReceiver<T>{
    buffer: VecDeque<T>,
    // buffersize: usize,
    rx: crossbeam::channel::Receiver<VecDeque<T>> // the inner sender has to handle collections of T!!
}

impl<T> BufferedReceiver<T>{
    pub fn new(rx: crossbeam::channel::Receiver<VecDeque<T>>, buffersize: usize)-> Self{
        let buffer = VecDeque::with_capacity(buffersize); // note: the buffer has variable size!can get larger than buffersize
        BufferedReceiver {buffer, rx}
    }

    pub fn recv(&mut self) -> Result<T, RecvError>{

        // if the buffer is empty, pull in more content
        if self.buffer.len() == 0 {
            // let newbuffer: VecDeque<T> = VecDeque::with_capacity(self.buffersize);
            println!("pulling in new buffer");
            self.buffer = self.rx.recv()?;
        }
        let item = self.buffer.pop_front().expect("cant happen, we chehck that the buffer is nonempty");
        Ok(item)
    }
}

pub fn buffered_unbounded<T>(buffersize: usize) -> (BufferedSender<T>, BufferedReceiver<T>){
    let (tx, rx) = unbounded();
    let buffered_tx = BufferedSender::new(tx, buffersize);
    let buffered_rx = BufferedReceiver::new(rx, buffersize);
    (buffered_tx, buffered_rx)
}

pub fn buffered_bounded<T>(bufersize: usize, bound:usize) -> (BufferedSender<T>, BufferedReceiver<T>){
    let (tx, rx) = bounded(bound);
    let buffersize = 1000;
    let buffered_tx = BufferedSender::new(tx, buffersize);
    let buffered_rx = BufferedReceiver::new(rx, buffersize);
    (buffered_tx, buffered_rx)
}

#[test]
fn test_buf_channel(){
    let (mut tx, mut rx) = buffered_unbounded(2);
    tx.send(1_usize).unwrap();
    tx.send(2_usize).unwrap();
    tx.send(3_usize).unwrap();
    println!("{:?}",tx.buffer);
    tx.flush().unwrap();

    let msg = rx.recv().unwrap();
    println!("recieved {:?}", msg);
    let msg = rx.recv().unwrap();
    println!("recieved {:?}", msg);    
    let msg = rx.recv().unwrap();
    println!("recieved {:?}", msg);   
}
